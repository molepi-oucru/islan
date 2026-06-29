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

### Preprocessing & Filtering (`pise pre-process`)
Run from the root directory:
- **Single-Sample Mode (using `uv`)**:
  ```bash
  uv run pise pre-process --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True --i5-mismatch 2 --min_len 20
  ```
- **Single-Sample Mode (using Conda)**:
  ```bash
  pise pre-process --config config/config.yaml /path/to/reads_1.fastq.gz /path/to/reads_2.fastq.gz --threads 8 --qc True --i5-mismatch 2 --min_len 20
  ```
- **Batch Mode for Multiple Samples (using `uv`)**:
  ```bash
  uv run pise pre-process --threads 8 --qc True /path/to/raw_reads_directory
  ```
- **Batch Mode for Multiple Samples (using Conda)**:
  ```bash
  pise pre-process --threads 8 --qc True /path/to/raw_reads_directory
  ```

### ASV Analysis (`pise asv-analysis`)
Run ASV analysis:
```bash
uv run pise asv-analysis --config config/config.yaml
```

### Pairing Reverse Reads (`pise pairing`)
Extract reverse reads matching your filtered forward reads or ID list:
- **Using forward FASTQ**:
  ```bash
  uv run pise pairing \
    -f results_asv/filtered_reads/sample_filtered_1.fastq.gz \
    -r /path/to/raw_reads_2.fastq.gz \
    -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
  ```
- **Using plain text list of IDs**:
  ```bash
  uv run pise pairing \
    -f /path/to/ids.txt \
    -r /path/to/raw_reads_2.fastq.gz \
    -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
  ```

### IS Mapping (`pise is-mapping`)
Map filtered reads to identify insertion sites on a reference genome.
- **WGS Mode**: Map raw forward/reverse reads:
  ```bash
  uv run pise is-mapping --mode wgs \
    -1 /path/to/sample_R1.fastq.gz -2 /path/to/sample_R2.fastq.gz \
    --is_query /path/to/IS.fasta \
    --ref_seq /path/to/reference.gbk \
    --output_dir results_wgs_mapping
  ```
- **Targeted Mode**: Map filtered and extracted HEAD/TAIL reads directly:
  ```bash
  uv run pise is-mapping --mode targeted \
    --filtered_forward results_asv/filtered_reads/sample_filtered_1.fastq.gz \
    --filtered_reverse results_asv/filtered_reads/sample_filtered_2.fastq.gz \
    --filtered_head results_asv/filtered_reads/sample_filtered_1_HEAD.fastq.gz \
    --filtered_tail results_asv/filtered_reads/sample_filtered_1_TAIL.fastq.gz \
    --ref_seq /path/to/reference.gbk \
    --min_endogenous_depth 50 \
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

