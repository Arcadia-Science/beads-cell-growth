"""Batch cell segmentation pipeline for ND2 microscopy images, deployed on Modal.

For each ND2 file, this script:
1. Loads all four channels (DIC, FITC, TRITC, DAPI).
2. Applies channel-specific preprocessing (background subtraction for fluorescence,
   percentile rescaling for all).
3. Segments cells from the DIC channel using a Cellpose model.
4. Saves per-cell morphological properties (CSV), a labeled segmentation mask (TIFF),
   and a multi-channel overlay with mask outlines (JPG) back to S3.

This script is designed to run remotely via Modal (https://modal.com) using GPU-backed
containers, and is NOT meant to be run locally. For a local example of the
segmentation workflow, see ``notebooks/segmentation-example.ipynb``.

Usage::

    modal run scripts/segment_cells.py
"""

from __future__ import annotations
from pathlib import Path

import modal

app = modal.App("segment-cells")

s3_mount = modal.CloudBucketMount(
    bucket_name="arcadia-microscopy-data",
    secret=modal.Secret.from_name("aws-s3-access-secret"),
)

dependencies = [
    "arcadia-microscopy-tools[all]>=0.2.5",
]
image = modal.Image.debian_slim().uv_pip_install(dependencies)


@app.function(image=image, volumes={"/s3": s3_mount})
def get_nd2_paths(prefix: str) -> list[str]:
    """Return sorted ND2 file paths under ``prefix`` on the mounted S3 bucket."""
    input_directory = Path("/s3") / prefix
    nd2_paths = sorted(input_directory.glob("*.nd2"))
    print(f"Input directory exists: {input_directory.exists()}")
    print(f"Found {len(nd2_paths)} ND2 files within {input_directory}")
    return [str(p) for p in nd2_paths]


@app.function(image=image, gpu="T4", volumes={"/s3": s3_mount}, max_containers=10)
def process_nd2_file(nd2_path_str: str, output_prefix: str) -> tuple[str, int]:
    """Segment cells in a single ND2 file and write results to S3.

    Preprocessing differs by channel type: fluorescence channels (FITC, TRITC,
    DAPI) receive difference-of-Gaussians background subtraction before
    percentile rescaling, while the DIC channel is only rescaled.

    Segmentation is performed on the preprocessed DIC channel because it
    provides the clearest cell morphology for Cellpose.  Raw (unprocessed)
    intensities are used when computing per-cell properties so that
    measurements reflect true fluorescence values.

    Outputs written to ``<output_prefix>/processed/`` on S3:
        - ``*_mask.tiff``       - 16-bit labeled segmentation mask
        - ``*_overlay.jpg``     - RGB overlay (DIC + fluorescence + mask outlines)
        - ``*_properties.csv``  - per-cell morphological / intensity properties

    Returns:
        Tuple of (input path, number of cells found). Returns 0 cells if
        segmentation finds nothing.
    """
    import shutil
    import tempfile
    import warnings

    import numpy as np
    import pandas as pd
    from arcadia_microscopy_tools import Channel, MicroscopyImage, overlay_channels
    from arcadia_microscopy_tools.channels import DAPI, DIC, FITC, TRITC
    from arcadia_microscopy_tools.masks import SegmentationMask
    from arcadia_microscopy_tools.model import SegmentationModel
    from arcadia_microscopy_tools.operations import rescale_by_percentile, subtract_background_dog
    from arcadia_microscopy_tools.pipeline import ImageOperation, Pipeline
    from arcadia_pycolor import HexCode
    from cellpose.utils import masks_to_outlines
    from skimage.exposure import rescale_intensity
    from skimage.io import imsave

    nd2_path = Path(nd2_path_str)
    channels = [DIC, FITC, TRITC, DAPI]
    image = MicroscopyImage.from_nd2_path(nd2_path, channels=channels)

    # Raw intensities are kept for computing per-cell properties later.
    intensity_image_dict = {
        channel: image.get_intensities_from_channel(channel) for channel in channels
    }

    # DIC only needs rescaling; fluorescence channels also need background
    # subtraction (difference-of-Gaussians) to remove uneven illumination.
    dic_pipeline = Pipeline(
        [ImageOperation(rescale_by_percentile, percentile_range=(0.1, 99.9))],
        preserve_dtype=False,
    )
    fluorescence_pipeline = Pipeline(
        [
            ImageOperation(subtract_background_dog, percentile=30),
            ImageOperation(rescale_by_percentile, percentile_range=(0.1, 99.9)),
        ],
        preserve_dtype=False,
    )
    pipelines = {
        DIC: dic_pipeline,
        FITC: fluorescence_pipeline,
        TRITC: fluorescence_pipeline,
        DAPI: fluorescence_pipeline,
    }

    model = SegmentationModel()

    preprocessed_intensities_dict = {}
    for channel in channels:
        pipeline = pipelines[channel]
        preprocessed = image.apply_pipeline(pipeline, channel).astype(float)
        preprocessed_intensities_dict[channel] = preprocessed

    try:
        cellpose_mask = model.segment(preprocessed_intensities_dict[DIC])
        segmentation_mask = SegmentationMask(
            cellpose_mask, intensity_image_dict=intensity_image_dict
        )
        if segmentation_mask.num_cells == 0:
            raise ValueError(f"No cells found for {nd2_path}")
    except ValueError:
        warnings.warn(f"No cells found for {nd2_path}", RuntimeWarning, stacklevel=2)
        return nd2_path_str, 0

    # Build an RGB overlay: DIC as grayscale background, fluorescence channels
    # in their assigned colors, and segmentation mask outlines in yellow.
    mask_channel = Channel("mask", color=HexCode("", "#ffff00"))
    mask_outlines = masks_to_outlines(segmentation_mask.label_image)
    background = preprocessed_intensities_dict[DIC]
    channel_intensities = {
        FITC: preprocessed_intensities_dict[FITC],
        TRITC: preprocessed_intensities_dict[TRITC],
        DAPI: preprocessed_intensities_dict[DAPI],
        mask_channel: mask_outlines,
    }
    overlay_rgb_float = overlay_channels(background, channel_intensities)
    overlay_rgb_8bit = rescale_intensity(overlay_rgb_float, out_range=np.ubyte)  # type: ignore

    pixel_size_um = image.metadata.image.channel_metadata_list[0].resolution.pixel_size_um
    cell_properties_um = segmentation_mask.convert_properties_to_microns(pixel_size_um)
    dataframe = pd.DataFrame(cell_properties_um)

    output_parent_path = Path(f"/s3/{output_prefix}/processed/")
    output_tiff_path = output_parent_path / f"{nd2_path.stem}_mask.tiff"
    output_jpg_path = output_parent_path / f"{nd2_path.stem}_overlay.jpg"
    output_csv_path = output_parent_path / f"{nd2_path.stem}_properties.csv"
    output_parent_path.mkdir(parents=True, exist_ok=True)

    # The S3 CloudBucketMount doesn't support direct writes from skimage/pandas,
    # so we write to a temp file first and then copy it to the mount path.
    segmentation_mask_16bit = segmentation_mask.label_image.astype(np.uint16)
    with tempfile.NamedTemporaryFile(suffix=".tiff", delete=False) as tmp:
        imsave(tmp.name, segmentation_mask_16bit, check_contrast=False)
        shutil.copy(tmp.name, str(output_tiff_path))
        Path(tmp.name).unlink()

    with tempfile.NamedTemporaryFile(suffix=".jpg", delete=False) as tmp:
        imsave(tmp.name, overlay_rgb_8bit)
        shutil.copy(tmp.name, str(output_jpg_path))
        Path(tmp.name).unlink()

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        dataframe.to_csv(tmp.name, index=False)
        shutil.copy(tmp.name, str(output_csv_path))
        Path(tmp.name).unlink()

    return nd2_path_str, segmentation_mask.num_cells


@app.local_entrypoint()
def main():
    """Local entrypoint invoked by ``modal run``.

    Change ``prefix`` to select which dataset to process. Each prefix
    corresponds to an S3 directory of ND2 files from a specific experiment
    (the comment after each prefix indicates the associated analysis script).
    """
    prefixes = [
        "Hina/Roman/2026-01-16/20260116_094944_372",  # 96_beads.py
        "Hina/Roman/2026-01-22/20260122_111821_521",  # ttubes_beads.py
        "Hina/Roman/2026-01-22/20260122_113404_129",  # 24_beads.py
        "Hina/Roman/2026-01-23/20260123_113447_096",  # supplements_beads.py
    ]

    prefix = prefixes[0]
    nd2_paths = get_nd2_paths.remote(prefix)

    print(f"Processing {len(nd2_paths)} files")
    for nd2_path, num_cells in process_nd2_file.starmap([(path, prefix) for path in nd2_paths]):
        nd2_filename = Path(nd2_path).name
        print(f"Segmented {num_cells} cells in {nd2_filename}")
