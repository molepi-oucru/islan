# Command Reference

This file compiles useful commands for running, testing, and managing the `pise` analysis pipeline.

## Installation

### Option A: Using `uv` (Recommended)
Global installation (registers commands to system PATH for execution from anywhere):
```bash
uv tool install --editable /path/to/PISE/
```

Local sync (creates local `.venv` environment):
```bash
uv sync --python 3.8
```

### Option B: Using Pip & Conda (Legacy)
Activate Conda environment:
```bash
conda activate py38_isseq
```
Install package in editable mode:
```bash
pip install -e .
```

## Running the Pipeline

### Orchestrator Pipeline
Run from the root directory:
- **Single-Sample Mode (using `uv`)**:
  ```bash
  uv run pise --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True
  ```
- **Single-Sample Mode (using Conda)**:
  ```bash
  pise --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True
  ```
- **Batch Mode for Multiple Samples (using `uv`)**:
  ```bash
  uv run pise --threads 8 --qc True /path/to/raw_reads_directory
  ```
- **Batch Mode for Multiple Samples (using Conda)**:
  ```bash
  pise --threads 8 --qc True /path/to/raw_reads_directory
  ```

### Deferred Reverse Reads Extraction
Extract reverse reads matching your filtered forward reads:
- **Using Python**:
  ```bash
  python pise/extract_pairs.py \
    -f results/filtered_reads/sample_index-TargetIS_1_full.fastq.gz \
    -r /path/to/raw_reads_2.fastq.gz \
    -o results/filtered_reads/sample_index-TargetIS_2_full.fastq.gz
  ```

## Running Tests
- **Using `uv`**:
  ```bash
  uv run python -m unittest discover -s tests
  ```
- **Using Conda**:
  ```bash
  python -m unittest discover -s tests
  ```
