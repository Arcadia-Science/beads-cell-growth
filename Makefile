MICROSCOPY_DIR ?= .

.PHONY: lint
lint:
	ruff check --exit-zero .

.PHONY: format
format:
	ruff check --fix .
	ruff format .

.PHONY: pre-commit
pre-commit:
	pre-commit run --all-files

.PHONY: compile-microscopy
compile-microscopy:
	python scripts/compile_microscopy.py -e 96-well     -d $(MICROSCOPY_DIR)/20260116_094944_372/processed
	python scripts/compile_microscopy.py -e ttubes       -d $(MICROSCOPY_DIR)/20260122_111821_521/processed
	python scripts/compile_microscopy.py -e 24-well      -d $(MICROSCOPY_DIR)/20260122_113404_129/processed
	python scripts/compile_microscopy.py -e supplements   -d $(MICROSCOPY_DIR)/20260123_113447_096/processed
