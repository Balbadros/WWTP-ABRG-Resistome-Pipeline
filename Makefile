PYTHON ?= python
CONFIG ?= config/example.yaml

setup:
	$(PYTHON) -m pip install -r requirements.txt

preprocess:
	$(PYTHON) -m wwtp_abrg.run --config $(CONFIG)

top30: preprocess

pca: preprocess

diff: preprocess

network: preprocess

figures: preprocess

all: preprocess
