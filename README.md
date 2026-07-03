# PISE: Processing IS-Seq Sequencing Data

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

### Option A: Using Conda & Mamba (Recommended - Fully Isolated)
This is the most convenient way to set up the pipeline, as it installs the Python environment, all Python libraries, and all compiled bioinformatic binary dependencies (`bwa`, `bwa-mem2`, `samtools`, `bedtools`, `blastn`, and `minimap2`) in a single command.

1. **Create and Activate Environment**:
   Run the following from the root of the repository to create the environment and register the `pise` tool in editable mode:
   ```bash
   mamba env create -f environment.yaml
   conda activate pise_env
   ```

2. **Verify Installation**:
   Verify that `pise` and all dependencies are correctly registered:
   ```bash
   pise --help
   samtools --version
   blastn -version
   minimap2 --version
   ```

### Option B: Using `uv` (Fast Python Setup)
If you already have the required command-line dependencies (`bwa`, `bwa-mem2`, `samtools`, `bedtools`, `blastn`, and `minimap2`) installed on your system `PATH`, you can use `uv` for a fast Python environment setup.

1. **Global Installation**:
   Install the pipeline globally as an editable tool using `uv`:
   ```bash
   uv tool install --editable .
   ```
   This registers the `pise` command directly in your user path.

2. **Local Project Environment**:
   Initialize a local virtual environment `.venv/` and sync dependencies:
   ```bash
   uv sync
   ```
   Then run commands prepended with `uv run`:
   ```bash
   uv run pise ...
   ```


## Running the Pipeline

The pipeline uses a modular subcommand-based command line interface:

### 1. Preprocessing & Filtering (`pise pre-process`)
Run the preprocessing step by passing your raw forward FASTQ reads file and the optional raw reverse FASTQ reads file:
```bash
pise pre-process \
  --config config/config.yaml \
  --threads 4 \
  --qc True \
  --i5-mismatch 2 \
  --min_len 20 \
  /path/to/raw_reads_1.fastq.gz \
  /path/to/raw_reads_2.fastq.gz
```
You can also run in **batch mode** by passing a directory of raw reads:
```bash
pise pre-process \
  --config config/config.yaml \
  --threads 4 \
  /path/to/raw_reads_directory
```

*Notes on Preprocessing:*
* QC & demultiplexing are controlled via `--qc True` (default) or `--qc False`. 
* Poly-N reads are dropped, and index demultiplexing is performed based on the 8-bp i5 index prefixes in `config/primers.fasta` with a mismatch tolerance specified by `--i5-mismatch` (default 2).
* Samples with expected index reads < 30% of total non-poly-N reads are classified as "Index failure" and automatically skipped in the downstream filtering.
* Valid samples are trimmed by 8 bp, length-filtered with `MIN_LEN_FORWARD`, and split into HEAD/TAIL files based on local alignment match.
* The processing statistics and classifications are saved to `pise_summary.tsv` in the output directory.

### 2. ASV Analysis (`pise asv-analysis`)
Runs Amplicon Sequence Variant analysis:
```bash
pise asv-analysis --config config/config.yaml
```

### 3. Extracting Reverse Reads (`pise pairing`)
Extract the reverse reads matching the filtered forward reads or a list of read IDs:
```bash
pise pairing \
  -f results_asv/filtered_reads/sample_filtered_1.fastq.gz \
  -r /path/to/raw_reads_2.fastq.gz \
  -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
```
Or with a plain text file containing one ID per line:
```bash
pise pairing \
  -f /path/to/id_list.txt \
  -r /path/to/raw_reads_2.fastq.gz \
  -o results_asv/filtered_reads/sample_filtered_2.fastq.gz
```

### 4. IS Mapping (`pise is-mapping`)
Map filtered reads to identify insertion sites on a reference genome.

- **WGS Mode**: Map raw forward/reverse reads:
  ```bash
  pise is-mapping \
    --reads /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz \
    --queries /path/to/IS.fasta \
    --reference /path/to/reference.gbk \
    --output_dir results_wgs_mapping
  ```
- **Targeted Mode**: Map filtered and extracted HEAD/TAIL reads directly:
  ```bash
  pise is-mapping --targeted \
    --head results_asv/filtered_reads/sample_filtered_1_HEAD.fastq.gz \
    --tail results_asv/filtered_reads/sample_filtered_1_TAIL.fastq.gz \
    --filtered_forward results_asv/filtered_reads/sample_filtered_1.fastq.gz \
    --filtered_reverse results_asv/filtered_reads/sample_filtered_2.fastq.gz \
    --reference /path/to/reference.gbk \
    --cutoff 6 \
    --min-mapq 30 \
    --flank-len 300 \
    --output_dir results_targeted_mapping
  ```

Outputs are sorted genome-wide and split into:
- `{sample}_table.tsv`: Main high-confidence insertion table (Known Pairs, Novel Pairs, and resolved Tandems) with complete flanking gene annotation details and median/IQR depth statistics.
- `{sample}_unpaired.tsv`: Unpaired singleton and off-target noise table with gene annotations omitted.

For details on the algorithm, see [Targeted_Mapping.md](docs/Targeted_Mapping.md).

## Running Tests
Run unit tests to verify package integrity:
```bash
uv run pytest
```
