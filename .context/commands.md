# Command Reference

This file compiles useful commands for running, testing, and managing the `islan` analysis pipeline.

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

### Preprocessing & Filtering (`islan pre-process`)
Run from the root directory:
- **Single-Sample Mode (using `uv`)**:
  ```bash
  uv run islan pre-process --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True --i5-mismatch 2 --min_len 20
  ```
- **Single-Sample Mode (using Conda)**:
  ```bash
  islan pre-process --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True --i5-mismatch 2 --min_len 20
  ```
- **Batch Mode for Multiple Samples (using `uv`)**:
  ```bash
  uv run islan pre-process --threads 8 --qc True /path/to/raw_reads_directory
  ```
- **Batch Mode for Multiple Samples (using Conda)**:
  ```bash
  islan pre-process --threads 8 --qc True /path/to/raw_reads_directory
  ```

### ASV Analysis (`islan asv-analysis`)
Run ASV analysis:
```bash
uv run islan asv-analysis --config config/config.yaml
```

### Pairing Reverse Reads (`islan pairing`)
Extract reverse reads matching your filtered forward reads or ID list:
- **Using forward FASTQ**:
  ```bash
  uv run islan pairing \
    -f results_asv/filtered_reads/sample_filtered_1.fastq.gz \
    -r /path/to/raw_reads_2.fastq.gz \
    -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
  ```
- **Using plain text list of IDs**:
  ```bash
  uv run islan pairing \
    -f /path/to/ids.txt \
    -r /path/to/raw_reads_2.fastq.gz \
    -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
  ```

### IS Mapping (`islan is-mapping`)
Map filtered reads to identify insertion sites on a reference genome.
- **WGS Mode**: Map raw forward/reverse reads:
  ```bash
  uv run islan is-mapping \
    --reads /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz \
    --queries /path/to/IS.fasta \
    --reference /path/to/reference.gbk \
    --output_dir results_wgs_mapping
  ```
- **Targeted Mode**: Map filtered and extracted HEAD/TAIL reads directly:
  ```bash
  uv run islan is-mapping --targeted \
    --head results_asv/filtered_reads/sample_filtered_1_HEAD.fastq.gz \
    --tail results_asv/filtered_reads/sample_filtered_1_TAIL.fastq.gz \
    --filtered_forward results_asv/filtered_reads/sample_filtered_1.fastq.gz \
    --filtered_reverse results_asv/filtered_reads/sample_filtered_2.fastq.gz \
    --reference /path/to/reference.gbk \
    --output_dir results_targeted_mapping
  ```


## Running Tests
- **Using `uv`**:
  ```bash
  uv run pytest
  ```
- **Using Conda**:
  ```bash
  pytest
  ```
