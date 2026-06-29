# Targeted IS-Seq Mapping Guide

This document describes the design, algorithms, parameter tuning, output formats, and mathematical logic behind the `pise is-mapping --targeted` pipeline.

---

## 1. Quick Start Example Command

To run the pipeline in targeted mode, supply the filtered forward reads, filtered reverse reads (ensuring adapters and index offsets are removed), and the extracted `HEAD`/`TAIL` read subsets:

```bash
uv run pise is-mapping --targeted \
  --head results_asv/filtered_reads/sample_1_HEAD.fastq.gz \
  --tail results_asv/filtered_reads/sample_1_TAIL.fastq.gz \
  --filtered_forward results_asv/filtered_reads/sample_filtered_1.fastq.gz \
  --filtered_reverse results_asv/filtered_reads/sample_filtered_2.fastq.gz \
  --reference references/genome.gb \
  --output_dir results_targeted \
  --cutoff 6 \
  --min-mapq 30 \
  --flank-len 300 \
  --threads 8
```

### Key CLI Parameters
*   `--targeted`: Activates targeted amplicon parsing mode, bypassing the standard WGS soft-clip parsing heuristics.
*   `--head` / `--tail`: Trimmed primer-specific read FASTQ files generated during `pise pre-process`.
*   `--filtered_forward` / `--filtered_reverse`: Trimmed paired-end reads (essential to prevent i5 indexes or adapters from affecting the soft-clipping and mapping coordinates).
*   `--cutoff`: Minimum read depth cutoff at each base position to consider it as part of a called flanking peak (default 6).
*   `--min-mapq`: Minimum mapping quality filter (default 30). This retains perfect multi-mappers (`MAPQ == 0`) and filters out weak cross-hybridizations (`0 < MAPQ < 30`).
*   `--flank-len`: Search window size around known reference copies (default 300 bp).

---

## 2. Why WGS Mode Should Not Be Used for Targeted IS-Seq

You must always run `pise is-mapping` with the `--targeted` flag when analyzing IS-Seq data. Running standard WGS mode (the original ISMapper algorithm) on targeted amplicon reads breaks down because of the difference in read structures.

### WGS shotgun vs. Amplicon Data Structure
- **WGS Shotgun Data**: DNA is fragmented randomly. Read pairs spanning an IS boundary consist of one read mapping entirely to the genomic sequence, and its mate spanning the boundary (partially mapping to the IS, and partially to the genome).
- **Targeted Amplicon Data**: Primers bind *inside* the IS element and point **outward** into the flanking genomic sequence. Thus:
  - **Read 1 (`_1`)** always starts with the IS primer, reading outward into the genomic flank.
  - **Read 2 (`_2`)** starts from the other end of the fragment, reading *backward* across the exact same genomic flank.

### The Failure of WGS Heuristics
The WGS algorithm (ISMapper) sorts shotgun reads into "Left Flank" and "Right Flank" buckets based on CIGAR soft-clipping and SAM orientation flags. When fed structured, outward-facing amplicon reads:
1. **Overlapping Flanks Artifact**: For a known/endogenous IS element, both reads in a pair sequence the exact same genomic region. Due to SAM flag sorting, the WGS algorithm accidentally sorts Read `_1` into the "Left Flank" bucket and Read `_2` into the "Right Flank" bucket.
2. **False Intersects & Data Loss**: The pipeline sees a "Left" and "Right" flank mapping to the same location and flags them as a false insertion. The WGS gap-size filter subsequently discards the hit because an overlap of ~1000+ bp violates novel insertion rules, resulting in silent data loss.

---

## 3. How Insertion Sites Are Detected

PISE resolves insertions using a **three-stage coordinate pairing logic** designed to identify known copies, characterize novel and tandem insertions, and isolate noise. 

### 3.1 Reference Mapping and Quality Filtering (MAPQ)
Before the pairing logic is executed, extracted flanking reads are aligned back to the reference sequence:
1. **Mapping with BWA**:
   - The flanking reads are mapped to the reference genome using `bwa mem`.
2. **Mapping Quality (MAPQ) Stream Filter**:
   - Alignments are filtered using an inline stream processor:
     ```bash
     awk '$5 == 0 || $5 >= min_mapq || $1 ~ /^@/'
     ```
   - **Why MAPQ == 0 is Retained (Multi-mappers)**: Perfect multi-mappers (alignments mapping to multiple locations with equal score, mapped with MAPQ=0 by BWA) are *deliberately retained*. This is crucial because the IS element may insert into multicopy genes (such as ribosomal RNA operons, homologous transposon regions, or duplicated genomic segments). Discarding MAPQ=0 reads would completely hide insertions in these repetitive regions.
   - **Filtering Low-Quality Alignments (0 < MAPQ < min_mapq)**: Reads with intermediate or low quality mapping (where BWA favors one position slightly but with low confidence) are discarded to avoid false positive calls arising from partial/spurious sequence homologies.
   - **High-Confidence Alignments (MAPQ >= min_mapq)**: Standard unique alignments (default threshold MAPQ >= 30) are preserved.
3. **Peak Calling and Coverage Depth Cutoff**:
   - The depth of mapped reads is calculated genome-wide.
   - To filter out low-coverage background noise, a minimal coverage filter is applied:
     ```bash
     awk '$4 >= {cutoff}' {left_cov} > {left_final_cov}
     ```
   - Only bases with a coverage depth greater than or equal to the `--cutoff` threshold (default 6) are preserved.
   - The remaining coordinates are merged using `bedtools merge -d {merging}` to call candidate flanking peaks. Any peak analyzed by the pairing logic is therefore supported by at least 6 reads.

### 3.2 Chromosome and Plasmid References
Bacteria often harbor plasmids in addition to the primary chromosome.
1. **Combined Reference Requirement**:
   - If an IS element is located on a plasmid, but only the chromosome is provided as the reference sequence, plasmid-derived reads may remain unmapped or mistakenly align with low-confidence to chromosomal regions sharing weak sequence homology.
   - To overcome this, a combined reference file (in a multi-fasta or multi-genbank file) containing both the chromosome and all plasmids can be used.
2. **Simultaneous Mapping (Parallel Competition)**:
   - When a combined reference containing both the chromosome and plasmids is provided, mapping is performed **simultaneously** rather than sequentially.
   - PISE parses all records from the input file, merges them into a single temporary reference FASTA, and builds a single combined BWA index.
   - During `bwa mem` execution, all reads are aligned against this combined index in a single run. This allows reads to compete across all replicons (chromosome and plasmids) at the same time. This is critical for evaluating correct mapping quality (MAPQ) and identifying true multi-mappers vs. unique plasmid/chromosomal insertions.
3. **Ambiguous Insertions (Similar/Duplicated Regions)**:
   - When identical sequence segments or homologous genes are present on both the chromosome and a plasmid, or duplicated on the chromosome:
     - The flanking reads will map to all copies with MAPQ=0.
     - The pipeline will report these as parallel hits on all homologous regions, reflecting the mathematical ambiguity of the insertion.
     - Users can resolve these by checking coverage depth (plasmids usually have higher copy numbers and thus higher coverage than the chromosome) or using long-read sequencing verification.

### 3.3 The Three-Stage Pairing Logic

#### Step 0: Chimera Filtering (Pre-filtering)
Before executing the three-stage pairing stages, the algorithm performs PCR chimera filtering on the raw peak pools:
- **Condition**: If a `HEAD` and `TAIL` peak overlap by 90%+ (`is_full_overlap`) and their coverage depth ratio exceeds `5.0`.
- **Action**: The minor peak is flagged as a chimera and immediately discarded from the active pools. It is excluded from all subsequent Stage 1, 2, or 3 pairing/singleton checks.

#### Stage 1: Resolve Known (Endogenous) IS Elements
1. **Target-Guided Scan**:
   - The reference genome is mapped against target elements carrying the `:FULL` suffix in the [targets.fasta](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/config/targets.fasta) file using `bwa mem -a`.
   - This scan determines the exact genomic coordinates and orientation of all known target copies present in the reference sequence.
2. **Peak Pairing**:
   - For each known IS element coordinate interval `[start, end]` and orientation:
     - In **Forward (`+`)** orientation: The algorithm searches for a `HEAD` peak near `start`, and a `TAIL` peak near `end` within a window of `--flank-len`.
     - In **Reverse (`-`)** orientation: The algorithm searches for a `TAIL` peak near `start`, and a `HEAD` peak near `end`.
   - If both peaks are detected, they are paired as a `Known Pair`, marked as `paired`, and removed from the active pools.
3. **Logging**: The pipeline reports the total known copies found on the reference and how many were successfully paired.

#### Stage 2: Resolve Novel Insertions and Tandems
Active un-paired peaks are processed using coordinate-based overlap check functions:
1. **Novel Same-Direction Tandem**:
   - Signature: Two insertions in the same direction (`++` or `--`) next to each other.
   - Coordinate check: `l_peak_1` (HEAD) partially overlaps with a central pair consisting of `r_peak_1` (TAIL) and `l_peak_2` (HEAD) that fully overlap, which in turn partially overlaps with `r_peak_2` (TAIL).
   - Reported as `Tandem Gap (same direction)`.
2. **Novel Opposite-Direction Tandem**:
   - Signature `+-` (pointing towards each other): Two `HEAD` peaks (`lp1`, `lp2`) both partially overlap with a single central `rpN` (TAIL). Reported as `Tandem Pair (+-)`.
   - Signature `-+` (pointing away from each other): Two `TAIL` peaks (`rp1`, `rp2`) both partially overlap with a single central `lpN` (HEAD). Reported as `Tandem Pair (-+)`.
3. **Novel Pair (with Target Site Duplications - TSD)**:
   - Signature: A remaining `HEAD` and `TAIL` peak partially overlap each other (overlap length usually $< 100$ bp).
   - Reported as `Novel Pair (TSD)`.
4. **Novel Pair (Standard)**:
   - Signature: A remaining `HEAD` and `TAIL` peak do not overlap but are within 100 bp of each other.
   - Reported as `Novel Pair`.

#### Stage 3: Resolve Singletons and Noise
1. **Off-Target Amplicon (Noise)**:
   - If a remaining `HEAD` and `TAIL` peak fully overlap (checked using `is_full_overlap`), they represent off-target PCR amplification.
   - Re-classified as `Off-Target Amplicon (Noise)` and moved to the unpaired TSV.
2. **Singletons**:
   - Remaining un-paired peaks are reported as `HEAD-only` or `TAIL-only` singletons in the unpaired TSV.

---

## 4. PCR Noise, Chimeras, and IS Orientation

### PCR Chimera Filtering (Depth Ratio check)
During PCR amplification of endogenous elements, high product concentrations can cause chimera artifacts (aborted extensions primer-matching a different flank).
*   **The Artifact**: A chimera maps to the exact same position as a true flank, creating a false `HEAD`/`TAIL` overlapping peak.
*   **The Solution**: True insertion overlaps (like TSDs) exhibit balanced depth. Chimeras show massive depth asymmetry.
*   **Rule**: If a `HEAD` and `TAIL` peak overlap by 90%+ (`is_full_overlap`), their depth ratio is calculated:
    $$\text{Ratio} = \frac{\max(\text{mean\_depth}_{\text{left}}, \text{mean\_depth}_{\text{right}})}{\min(\text{mean\_depth}_{\text{left}}, \text{mean\_depth}_{\text{right}})}$$
*   If the ratio is $> 5.0$, the minor peak is flagged as a chimera and discarded.

### Resolving IS Orientation
Based on outward-facing primers, the orientation of a paired insertion is deduced using their coordinates:
*   **Forward (`F` / `+`) Orientation**: The 5' end (`HEAD`) is on the genomic left, and the 3' end (`TAIL`) is on the genomic right (`HEAD_coordinate < TAIL_coordinate`).
*   **Reverse (`R` / `-`) Orientation**: The 3' end (`TAIL`) is on the genomic left, and the 5' end (`HEAD`) is on the genomic right (`TAIL_coordinate < HEAD_coordinate`).

---

## 5. Output Data Format & Interpretation

PISE outputs two main TSV files:
1. **`{sample}_table.tsv`**: Paired insertions (Known, Novel, and resolved Tandems) with flanking gene annotations.
2. **`{sample}_unpaired.tsv`**: Unpaired singletons and off-target noise (with gene annotations omitted).

### Output Column Descriptions

| Column Name | Description |
| :--- | :--- |
| `region` | Unique index sequence (`region_1`, `region_2` etc. in table; `unpaired_1`, `unpaired_2` in unpaired). |
| `contig` | Name of the reference sequence (e.g. chromosome or plasmid) where the peaks mapped. |
| `orientation` | Predicted orientation of the insertion: `F` (Forward), `R` (Reverse), or `?` (for singletons/noise). |
| `x` | Left-most boundary of the insertion site. |
| `y` | Right-most boundary of the insertion site. |
| `gap` | Gap distance between insertion site boundaries. Positive for non-overlapping gaps; negative for overlaps (TSDs). |
| `call` | Classification of the hit: `known`, `novel`, `novel (TSD)`. |
| `left_pos` | Chromosome coordinate range of the left flanking region. |
| `right_pos` | Chromosome coordinate range of the right flanking region. |
| `left_depth_median` | Median coverage depth across the left flanking peak. |
| `left_depth_iqr` | Interquartile Range (first and third quartiles `Q1-Q3`) of the left flanking peak's depth. |
| `right_depth_median` | Median coverage depth across the right flanking peak. |
| `right_depth_iqr` | Interquartile Range (`Q1-Q3`) of the right flanking peak's depth. |
| `left_gene` | Locus tag of the closest gene upstream of the insertion. |
| `left_description` | Product/description of the upstream gene. |
| `left_strand` | Coding strand of the upstream gene (`1` or `-1`). |
| `left_distance` | Distance (bp) from the insertion boundary `x` to the upstream gene. Reported as `0` if intragenic (interrupted). |
| `right_gene` | Locus tag of the closest gene downstream of the insertion. |
| `right_description` | Product/description of the downstream gene. |
| `right_strand` | Coding strand of the downstream gene (`1` or `-1`). |
| `right_distance` | Distance (bp) from the insertion boundary `y` to the downstream gene. |
| `gene_interruption`| Evaluates as `True` if the insertion is located within the coding sequence (CDS) of either flanking gene. |
