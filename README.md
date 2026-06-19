# pise: IS-Seq Sequencing Data Analysis Pipeline

A modular Python package and command-line pipeline for preprocessing and analyzing Insertion Sequence Sequencing (IS-Seq) reads.

## Repository Layout
```
PISE/
├── .context/               # AI Agent context directory
│   ├── algorithms.md
│   ├── architecture.md
│   ├── commands.md
│   ├── README.md
│   └── ...
├── config/                 # Configurations templates
│   ├── config.yaml         # Configuration file
│   └── primers.fasta       # Primers database
├── pise/                   # Core Python package
│   ├── __init__.py
│   ├── main.py             # Main orchestrator entry point
│   ├── qc.py               # Quality Control & demultiplexing logic
│   ├── filter_reads.py     # Parallelized read filtering algorithm
│   ├── extract_pairs.py    # Deferred matching reverse reads extraction
│   └── split_fastq.py      # Utility to split FASTQ files
├── tests/                  # Automated tests
│   └── test_preprocess.py  # Preprocessing unit test
├── environment.yaml        # Conda environment definition
├── pyproject.toml          # PEP 517 build configuration
├── setup.py                # Legacy setup compatibility file
└── README.md               # User manual (this file)
```

## Setup & Installation

### Option A: Using `uv` (Recommended - Blazingly Fast)
`uv` is an extremely fast Python package manager that handles virtual environment creation and package installation automatically.

1. **Global Installation (Run from anywhere)**:
   You can install the pipeline as a global CLI tool using `uv`. This automatically places the `pise` executable on your system's `PATH`:
   ```bash
   uv tool install --editable /path/to/PISE/
   ```
   Once installed, you can call `pise` directly from any directory in your terminal!

2. **Alternative: Local Project Environment**:
   Initialize a local virtual environment `.venv/` and install dependencies in editable mode:
   ```bash
   uv sync --python 3.8
   ```
   Then, prepend commands with `uv run` to execute them:
   ```bash
   uv run pise ...
   ```

### Option B: Using Conda & Pip (Legacy)

1. **Environment Setup**:
   Create and activate your Python/Conda environment (requires Python >= 3.8 and Biopython):
   ```bash
   conda env create -f environment.yaml
   conda activate py38_isseq
   ```

2. **Package Installation**:
   Install the package in editable mode from the repository root:
   ```bash
   pip install -e .
   ```
This registers the CLI command `pise` directly in your environment.

## Running the Pipeline

The pipeline processes forward reads first (QC -> Demultiplexing -> Filtering) and defers reverse read matching.

### 1. Preprocessing & Filtering
Run the orchestrator by passing your raw forward FASTQ reads file and the optional raw reverse FASTQ reads file:
```bash
pise \
  --config config/config.yaml \
  --threads 4 \
  --qc True \
  /path/to/raw_reads_1.fastq.gz \
  /path/to/raw_reads_2.fastq.gz
```
You can also run in **batch mode** by passing a directory of raw reads:
```bash
pise \
  --config config/config.yaml \
  --threads 4 \
  /path/to/raw_reads_directory
```

*Note on QC / Demultiplexing:*
* QC is controlled via `--qc True` (default) or `--qc False`. 
* Poly-N reads are dropped, and index demultiplexing is performed based on the 8-bp i5 index prefixes in `config/primers.fasta`.
* The filtering step will automatically skip files containing unknown indices (`index-unknown`).

### 2. Extracting Reverse Reads
After filtering, you can extract the reverse reads matching the filtered forward reads using `extract_pairs.py`:
```bash
python pise/extract_pairs.py \
  -f results/filtered_reads/sample_index-TargetIS_1_full.fastq.gz \
  -r /path/to/raw_reads_2.fastq.gz \
  -o results/filtered_reads/sample_index-TargetIS_2_full.fastq.gz
```

## Running Tests
Run unit tests to verify package integrity:
```bash
python -m unittest discover -s tests
```
