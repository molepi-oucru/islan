# IS-Seq Pipeline Context Overview

Welcome! This is the context directory designed for Gemini (AI Agent) to understand the architecture, purpose, and coding guidelines of the IS-Seq preprocessing and analysis pipeline (`pise`).

## Project Purpose
This repository implements `pise`, a modular Python package and pipeline for analyzing paired-end IS-Seq (Insertion Sequence Sequencing) data. IS-Seq sequences the junctions between insertion sequence (IS) elements and the surrounding bacterial/host genomes to discover transposition events and quantify their abundance.

## Modules Overview
1. **Quality Control & Demultiplexing (`pise/qc.py`)**: Performs fast 4-line parsing to count reads, drop poly-N reads, and demultiplex forward reads based on the 8-bp i5 index prefixes in `config/primers.fasta`. Writes demultiplexed reads to the `demux_reads/` subfolder.
2. **Reads Filter (`pise/filter_reads.py`)**: Implements the parallelized reads filtering algorithm using `multiprocessing.Pool` (local alignments using `PairwiseAligner`). It trims the 8-bp index, performs length filtering using `MIN_LEN_FORWARD`, and extracts the exact `len(primer)` matching prefix for HEAD/TAIL outputs. Writes filtered/HEAD/TAIL outputs to the `filtered_reads/` subfolder.
3. **Reverse Reads Pairing (`pise/extract_pairs.py`)**: A decoupled extraction script that extracts raw reverse reads matching filtered forward reads or a list of read IDs.
4. **IS Mapping (`pise/ismapper_adapted/`)**: Maps preprocessed reads to reference genome using BWA mem to identify insertion sites.
   - **WGS Mode**: Identifies insertion sites from shotgun reads using soft-clip heuristics.
   - **Targeted Mode**: Custom coordinate-based pairing logic utilizing 3 stages (known target scan, novel/tandem overlaps pairing, noise/singleton isolation).
   
For detailed descriptions of the underlying matching criteria, algorithms, and parallelization, see [algorithms.md](file:///.context/algorithms.md).

## Pipeline Orchestrator
- **Orchestrator (`pise/main.py` / `pise` command)**: A unified Python entrypoint that parses configuration templates and routes subcommands.
  - **`pre-process` Subcommand**: Coordinated by `pise/preprocess.py`, runs parallel QC over samples, dynamically determines index status (OK/warning/failure), skips `Index failure` samples, and performs read filtering/trimming.
  - **`pairing` Subcommand**: Delegates to `pise/extract_pairs.py` to extract matching reverse reads.
  - **`is-mapping` Subcommand**: Custom junction calling and insertion site characterization utilizing BWA mem, bedtools, and Biopython GenBank parsing.
  - **`asv-analysis` Subcommand**: A placeholder for downstream amplicon sequence variant analysis.
