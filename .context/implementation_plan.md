# Implementation Plan - Preprocessing Split & Filtering Parallelization

We want to:
1. Split the filtering step logic from `pise/preprocess_reads.py` to a new module `pise/filter_reads.py`. This leaves `preprocess_reads.py` as the overall coordinator for the preprocessing stage, ready for future optional QC and Trimming steps.
2. Parallelize the slow forward reads primer-filtering step in `pise/filter_reads.py` using Python's `multiprocessing` library to optimize processing speed on multi-core systems (e.g., the 10-core environment).

---

## User Review Required

> [!IMPORTANT]
> - **Execution Order Preservation**: We will use `multiprocessing.Pool.imap` to process batches of reads. This ensures that the output FASTQ files maintain the exact same read order as the input files, which is critical for downstream paired-end alignment tools.
> - **Thread Configuration**: We will expose a `threads` option in `config/config.yaml` under the `preprocessing` block, and also add a command-line argument `-t`/`--threads` to both `main.py` and `preprocess_reads.py`.
> - **Zero-dependency Multiprocessing**: This optimization uses standard Python libraries (`multiprocessing`, `functools`) and requires no additional packages.

---

## Proposed Changes

### Configuration Layer

#### [MODIFY] [config.yaml](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/config/config.yaml)
- Add `threads` config option under the `preprocessing` block:
  ```yaml
  preprocessing:
    min_cov: 0.98
    min_identity: 0.9
    min_len: 20
    index_len: 8
    len_actual_primer: null
    threads: 4 # Number of parallel processes to use for read filtering
  ```

---

### Pipeline Code Layer

#### [NEW] [filter_reads.py](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/pise/filter_reads.py)
- Create a dedicated file for the read filtering step.
- Implement `process_batch_worker` at the module level. This function:
  - Initializes a local `PairwiseAligner` inside the worker process (avoiding sharing/serialization issues).
  - Processes a batch of read tuples `(record_id, sequence, phred_quality, description)`.
  - Performs local alignment and checks full/partial criteria.
  - Collects local batch statistics and returns filtered read tuples.
- Refactor `process_forward_reads` to:
  - Add a `threads` argument.
  - If `threads <= 1`, run sequentially in the main thread (avoiding multiprocessing overhead).
  - If `threads > 1`, chunk the input reads into batches of 5000 and submit them to a `multiprocessing.Pool` via `imap` (preserving input order).
  - Reconstruct `Bio.SeqRecord` objects and write them to output files in the main process sequentially as batches complete.
- Keep the `process_reverse_reads`, `load_primers`, and `main()` parser logic inside `filter_reads.py`.

#### [MODIFY] [preprocess_reads.py](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/pise/preprocess_reads.py)
- Refactor this file to act as the overall coordinator.
- Import `load_primers`, `process_forward_reads`, and `process_reverse_reads` from `pise.filter_reads` to expose them so that `main.py` and any other existing callers do not break.
- Update the CLI `main()` entrypoint to call `filter_reads.main()`, or call `filter_reads` logic with additional parameters.

#### [MODIFY] [main.py](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/pise/main.py)
- Read `threads` from the `preprocessing` config block (defaulting to 4).
- Add a `-t`/`--threads` CLI argument. If provided, override the config-defined threads.
- Pass the determined `threads` to the preprocessing/filtering calls.

---

## Verification Plan

### Automated Tests
- Run `conda run -n py38 python -m unittest discover -s tests` to ensure imports and basic functionality work without errors.

### Manual Verification
- **Correctness Check**:
  Run the parallelized pipeline using the existing test data:
  ```bash
  conda run -n py38 python pise/main.py /data/SiNguyen/1.SIXTEEN/IS-SEQ/miseq_20250518/raw/278-IS1R_1.fastq.gz /data/SiNguyen/1.SIXTEEN/IS-SEQ/miseq_20250518/raw/278-IS1R_2.fastq.gz
  ```
  Compare the newly generated TSV summary files and gzipped FASTQ file contents against the backed-up originals using `diff` or md5 checksums to verify that they are mathematically identical.
- **Performance Measurement**:
  Measure the execution time of the preprocessing stage with `threads: 1` vs `threads: 8` to quantify the speedup.
