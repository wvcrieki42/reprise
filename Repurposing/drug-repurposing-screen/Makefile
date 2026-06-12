.PHONY: install demo test full clean

install:
	pip install -r requirements.txt

# Run the screen on the bundled sample data (no downloads needed)
demo:
	PYTHONPATH=src python -m repurpose.cli --config config.yaml --top 25

# End-to-end smoke test
test:
	PYTHONPATH=src python tests/test_smoke.py

# Run on the full downloaded databases (edit config.yaml: mode -> full, repoint paths)
full:
	PYTHONPATH=src python -m repurpose.cli --config config.full.yaml --top 50

clean:
	rm -rf output .cache **/__pycache__ src/**/__pycache__
