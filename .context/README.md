# IS-Seq Pipeline Context Overview

Welcome! This is the context directory designed for Gemini (AI Agent) to understand the architecture, purpose, and coding guidelines of the IS-Seq preprocessing and analysis pipeline (`pise`).

## Project Purpose
This repository implements `pise`, a modular Python package and pipeline for analyzing paired-end IS-Seq (Insertion Sequence Sequencing) data. IS-Seq sequences the junctions between insertion sequence (IS) elements and the surrounding bacterial/host genomes to discover transposition events and quantify their abundance.

## Modules Overview
1. **Quality Control & Demultiplexing (`pise/qc.py`)**: Performs fast 4-line parsing to count reads, drop poly-N reads, and demultiplex forward reads based on the 8-bp i5 index prefixes in `config/primers.fasta`.
2. **Reads Filter (`pise/filter_reads.py`)**: Implements the parallelized reads filtering algorithm using `multiprocessing.Pool` (local local alignments using `PairwiseAligner`). It classifies reads into full, partial, and binned outputs.
3. **Reverse Reads Pairing (`pise/extract_pairs.py`)**: A decoupled extraction script that extracts raw reverse reads matching filtered forward reads.
4. **Alignment (Future)**: Maps preprocessed reads to a reference genome using `bwa-mem2`.
5. **Junction Calling (Future)**: Analyzes mapped reads to identify insertion coordinates and calculate site abundance.

For detailed descriptions of the underlying matching criteria, algorithms, and parallelization, see [algorithms.md](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/.context/algorithms.md).

## Pipeline Orchestrator
- **Orchestrator (`pise/main.py` / `pise` command)**: A unified Python entrypoint that parses configuration templates. It coordinates execution in two phases:
  - **Phase 1: Parallel QC**: Spawns multiple samples in parallel using a `ProcessPoolExecutor` to run QC and demultiplexing.
  - **Phase 2: Filtering**: Runs filtering sequentially per sample but parallelizes execution internally using a `multiprocessing.Pool`. Automatically skips filtering for unknown index files.
