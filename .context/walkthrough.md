# Walkthrough - Modular Decoupled Preprocessing & Sample-Level Parallel QC

We have successfully restructured the reads preprocessing pipeline into a highly modular, decoupled, and parallelized architecture.

## Architecture Highlights

1. **Integrated QC & Demultiplexing (`pise/qc.py`)**
   - Implements fast 4-line parsing for forward reads.
   - Discards poly-N reads first.
   - Demultiplexes reads into up to 5 bins based on the 8-bp i5 index sequences defined in `config/primers.fasta`.
   - Generates interactive Plotly-based HTML quality control reports.

2. **Parallelized Reads Filtering (`pise/filter_reads.py`)**
   - Aligns reads against HEAD/TAIL primers using `Bio.Align.PairwiseAligner`.
   - Groups reads into batches of 5000 and distributes them to a `multiprocessing.Pool` of worker processes.
   - Preserves read order in output files by using order-guaranteed `.imap`.
   - Bypasses alignment/filtering entirely for reads demultiplexed as `index-unknown` to maximize CPU efficiency.

3. **Deferred Reverse Read Pairing (`pise/extract_pairs.py`)**
   - Decoupled reverse read processing from the main filtering loop.
   - Matches and extracts reverse reads matching the filtered forward reads in a single fast, memory-efficient pass after forward filtering is complete.

4. **Two-Stage Orchestrator (`pise/main.py`)**
   - **Phase 1 (Parallel QC)**: Spawns parallel QC processes across samples using `concurrent.futures.ProcessPoolExecutor` (taking advantage of multiple CPU cores for raw I/O/compression).
   - **Phase 2 (Internally Parallel Filtering)**: Filters demultiplexed files sequentially per sample, with each file utilizing a `multiprocessing.Pool` of size `threads`. This layout avoids nested process pools and locks.

5. **Clean Command-Line Interface**
   - Cleaned up the `--qc` CLI option and `qc:` config option from a string choice to a simpler boolean switch (`True`/`False`), keeping backwards compatibility with old yaml choices.
   - Relocated helper functions out of function scopes (avoiding nested definitions) for cleaner coding patterns.

## Verification

- All unit tests in `tests/test_preprocess.py` have been updated and pass successfully.
- End-to-end runs complete successfully with all expected outputs: HTML reports, demultiplexed files, filtered forward files, extracted reverse files, and summary files.
