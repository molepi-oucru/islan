# System Instructions for Gemini

When working on this repository, please adhere strictly to the following instructions:

## Code Quality and Environment
1. **Target Environment**: Python 3.8+ (using `uv` package manager/`.venv` virtual environment or Conda environment `py38`).
2. **Coding Style**: Strictly follow PEP 8 standards. Use meaningful variable names, docstrings, and comments.
3. **Libraries**:
   - For sequence operations, use **Biopython** (`Bio.Seq`, `Bio.SeqIO`, `Bio.Align`).
   - For file manipulation, use standard library (`os`, `sys`, `gzip`, `argparse`, `functools`).
   - For invoking external bioinformatics tools (like `bwa-mem2`), use Python's `subprocess` with `check_call` or `run` with safe argument lists (no shell=True unless necessary).
   - For parallel processing, utilize Python's standard `multiprocessing` library with batch-based chunking to minimize IPC overhead and maintain read pairing order.

## Bioinformatics Best Practices
1. **Gzipped Files**: Always read/write Fastq files directly in gzipped format using `gzip.open(..., "rt")` or `gzip.open(..., "wt")`. Do not unzip files on disk.
2. **Sequence Alignment**:
   - Use `Bio.Align.PairwiseAligner` in `'local'` mode for primer alignment.
   - Use `-1.0` penalty for mismatches and gap opens/extensions to avoid spuriously matching long unrelated regions.
3. **Error Handling**: Gracefully handle missing files, corrupted FASTQ formats, and missing configuration keys, printing informative messages to `sys.stderr` and exiting with non-zero exit codes.
