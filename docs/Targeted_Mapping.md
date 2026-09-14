# Targeted IS-Seq Mapping Guide

This document describes the design, algorithms, parameter tuning, output formats, sequence compositions, and mathematical logic behind the `islan is-mapping --targeted` pipeline.

---

## 1. Quick Start Example Command

To run the pipeline in targeted mode, supply the filtered forward reads, filtered reverse reads (ensuring adapters and index offsets are removed), and the extracted `HEAD`/`TAIL` read subsets:

```bash
uv run islan is-mapping --targeted \
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

Alternatively, to run the pipeline with forward reads only (single-end analysis), add the `--forward-only` flag (in this mode, `--filtered_reverse` is not required):

```bash
uv run islan is-mapping --targeted \
  --head results_asv/filtered_reads/sample_1_HEAD.fastq.gz \
  --tail results_asv/filtered_reads/sample_1_TAIL.fastq.gz \
  --filtered_forward results_asv/filtered_reads/sample_filtered_1.fastq.gz \
  --forward-only \
  --reference references/genome.gb \
  --output_dir results_targeted_forward_only \
  --cutoff 6 \
  --min-mapq 30 \
  --flank-len 300 \
  --threads 8
```

### Key CLI Parameters
*   `--targeted`: Activates targeted amplicon parsing mode, bypassing the standard WGS soft-clip parsing heuristics.
*   `--head`: Path to `{sample}_filtered_1_HEAD.fastq.gz`, containing reads matching the 5' IS terminus sequence (`HEAD`).
*   `--tail`: Path to `{sample}_filtered_1_TAIL.fastq.gz`, containing reads matching the 3' IS terminus sequence (`TAIL`).
*   `--filtered_forward`: Path to `{sample}_filtered_1.fastq.gz` (i5 index-trimmed Read 1 files).
*   `--filtered_reverse`: Path to `{sample}_filtered_2.fastq.gz` (paired Read 2 files).
*   `--forward-only`: Runs the mapping and analysis with forward reads (Read 1) only (single-end). In targeted mode, this ignores/skips reverse reads and doesn't require `--filtered_reverse`. In WGS mode, this configures single-end BWA mapping.
*   `--cutoff`: Minimum read depth cutoff at each base position to consider it as part of a called flanking peak (default 6).
*   `--min-mapq`: Minimum mapping quality filter (default 30). This retains perfect multi-mappers (`MAPQ == 0`) and filters out weak cross-hybridizations (`0 < MAPQ < 30`).
*   `--flank-len`: Search window size around known reference copies (default 300 bp).

---

## 2. Terminology & Sequence Composition of Input/Output Files

### 2.1 Distinguishing PCR Primers (`P_UP` / `P_DOWN`) vs. IS Termini (`HEAD` / `TAIL`)

To avoid confusion, ISLAN makes a strict distinction between the physical PCR primers and the IS element terminal regions:

*   **Physical PCR Primers (`P_UP` & `P_DOWN`)**:
    *   `P_UP` and `P_DOWN` are the actual physical PCR primers used in IS-Seq library preparation.
    *   Each primer sequence in the target definition begins with an 8 bp i5 index prefix followed by the actual primer annealing sequence of length `len_actual_primer` (~20–30 bp).
    *   `P_UP` anneals near the 5' end of the IS element, pointing **outward** into the upstream genomic flank.
    *   `P_DOWN` anneals near the 3' end of the IS element, pointing **outward** into the downstream genomic flank.

*   **IS Element Termini (`HEAD` & `TAIL`)**:
    *   `HEAD` and `TAIL` refer to the **structural 5' and 3' terminal regions** of the IS element body (typically ~80–100 bp sequences defined in `targets.fasta`).
    *   **`HEAD` (5' IS Terminus)**: The 5' structural end of the IS element. The sequence of a forward read matching `HEAD` starts at position 0 with the `P_UP` primer sequence (~`len_actual_primer` bp) and extends through the rest of the ~80–100 bp 5' IS terminal sequence before entering the host genomic sequence.
        *   During reference mapping in forward (`F` / `+`) orientation, `HEAD`-derived flanking reads map to the **left (upstream)** genomic flank of the insertion site.
    *   **`TAIL` (3' IS Terminus)**: The 3' structural end of the IS element. The sequence of a forward read matching `TAIL` starts at position 0 with the `P_DOWN` primer sequence (~`len_actual_primer` bp) and extends through the rest of the ~80–100 bp 3' IS terminal sequence before entering the host genomic sequence.
        *   During reference mapping in forward (`F` / `+`) orientation, `TAIL`-derived flanking reads map to the **right (downstream)** genomic flank of the insertion site.

---

### 2.2 Detailed Sequence Composition of Filtered Fastq Files

During `islan pre-process` (or `islan filter`), raw MiSeq FASTQ files are demultiplexed, trimmed, and segregated into specific intermediate files:

| File Name | Read Type | Trimming / Filtering Applied | Exact Sequence Composition |
| :--- | :--- | :--- | :--- |
| **`{sample}_filtered_1.fastq.gz`** | Forward (Read 1) | 1. Min length filter (`min_len_forward`).<br>2. Optional i5 index mismatch check.<br>3. 5' i5 index trimmed (`8 bp`). | **[IS Terminus Sequence (`HEAD` or `TAIL` starting with `P_UP`/`P_DOWN`)] + [Genomic Flanking Sequence]**<br>The read begins at position 0 with the IS terminal sequence (starting with the ~20–30 bp `P_UP` or `P_DOWN` primer sequence), extending through the IS terminal boundary into the host chromosome flank. |
| **`{sample}_filtered_2.fastq.gz`** | Reverse (Read 2) | Extracted during `islan pairing` by matching unique read IDs from `filtered_1`. | **[Reverse Genomic Flanking Sequence]**<br>Sequences from the opposite end of the PCR fragment, reading backward across the exact same genomic flank towards the IS element. |
| **`{sample}_filtered_1_HEAD.fastq.gz`** | Forward (Read 1) | Classified by local alignment matching the 5' IS terminus (`HEAD`). | **[HEAD IS Terminus Sequence Only]**<br>The sequence in each FASTQ record is truncated to **only the HEAD terminal sequence itself** (length = `head_primer_len`, ~80–100 bp). This file serves as a lookup registry of Read IDs sequencing the 5' (`HEAD`) IS boundary. |
| **`{sample}_filtered_1_TAIL.fastq.gz`** | Forward (Read 1) | Classified by local alignment matching the 3' IS terminus (`TAIL`). | **[TAIL IS Terminus Sequence Only]**<br>The sequence in each FASTQ record is truncated to **only the TAIL terminal sequence itself** (length = `tail_primer_len`, ~80–100 bp). This file serves as a lookup registry of Read IDs sequencing the 3' (`TAIL`) IS boundary. |

---

### 2.3 How Flanking Sequences Are Prepared for Reference Mapping

When `islan is-mapping --targeted` is executed, the pipeline processes the files as follows:

1. **ID & Terminal Length Extraction**: Read IDs and exact IS terminal sequence lengths (`head_primer_len` and `tail_primer_len`) are extracted from `{sample}_filtered_1_HEAD.fastq.gz` and `{sample}_filtered_1_TAIL.fastq.gz`.
2. **IS Sequence Trimming for Pure Genomic Flanks**:
   - Reads matching `HEAD` IDs are fetched from `{sample}_filtered_1.fastq.gz` and trimmed by `head_primer_len` (`seq[head_primer_len:]`). This strips away the entire IS `HEAD` terminal sequence (including `P_UP`), leaving **pure 5' genomic flanking sequence**, which is written to the temporary `left_final.fastq` pool.
   - Reads matching `TAIL` IDs are fetched from `{sample}_filtered_1.fastq.gz` and trimmed by `tail_primer_len` (`seq[tail_primer_len:]`). This strips away the entire IS `TAIL` terminal sequence (including `P_DOWN`), leaving **pure 3' genomic flanking sequence**, which is written to the temporary `right_final.fastq` pool.
3. **Paired Reverse Integration (unless `--forward-only`)**:
   - Reverse reads (Read 2) corresponding to `HEAD` IDs are concatenated into `left_final.fastq`.
   - Reverse reads (Read 2) corresponding to `TAIL` IDs are concatenated into `right_final.fastq`.
4. **Alignment**: `left_final.fastq` (5' genomic flanks) and `right_final.fastq` (3' genomic flanks) are mapped independently to the reference genome using `bwa mem`.

---

## 3. Why WGS Mode Should Not Be Used for Targeted IS-Seq

You must always run `islan is-mapping` with the `--targeted` flag when analyzing IS-Seq data. Running standard WGS mode (the original ISMapper algorithm) on targeted amplicon reads breaks down because of the difference in read structures.

### WGS shotgun vs. Amplicon Data Structure
- **WGS Shotgun Data**: DNA is fragmented randomly. Read pairs spanning an IS boundary consist of one read mapping entirely to the genomic sequence, and its mate spanning the boundary (partially mapping to the IS, and partially to the genome).
- **Targeted Amplicon Data**: Primers (`P_UP` / `P_DOWN`) bind *inside* the IS element termini (`HEAD` / `TAIL`) and point **outward** into the flanking genomic sequence. Thus:
  - **Read 1 (`_1`)** always starts with the IS terminal sequence (`HEAD` or `TAIL` starting with `P_UP`/`P_DOWN`), reading outward into the genomic flank.
  - **Read 2 (`_2`)** starts from the other end of the fragment, reading *backward* across the exact same genomic flank.

### The Failure of WGS Heuristics
The WGS algorithm (ISMapper) sorts shotgun reads into "Left Flank" and "Right Flank" buckets based on CIGAR soft-clipping and SAM orientation flags. When fed structured, outward-facing amplicon reads:
1. **Overlapping Flanks Artifact**: For a known/endogenous IS element, both reads in a pair sequence the exact same genomic region. Due to SAM flag sorting, the WGS algorithm accidentally sorts Read `_1` into the "Left Flank" bucket and Read `_2` into the "Right Flank" bucket.
2. **False Intersects & Data Loss**: The pipeline sees a "Left" and "Right" flank mapping to the same location and flags them as a false insertion. The WGS gap-size filter subsequently discards the hit because an overlap of ~1000+ bp violates novel insertion rules, resulting in silent data loss.

---

## 4. How Insertion Sites Are Detected

ISLAN resolves insertions using a **three-stage coordinate pairing logic** designed to identify known copies, characterize novel and tandem insertions, and isolate noise. 

### 4.1 Reference Mapping and Quality Filtering (MAPQ)
Before the pairing logic is executed, extracted flanking reads are aligned back to the reference sequence:
1. **Mapping with BWA**:
   - The flanking reads (`left_final.fastq` containing 5'/HEAD genomic flanks and `right_final.fastq` containing 3'/TAIL genomic flanks) are mapped to the reference genome using `bwa mem`.
2. **Mapping Quality (MAPQ) Stream Filter**:
   - Alignments are filtered using an inline stream processor:
     ```bash
     awk '$5 == 0 || $5 >= min_mapq || $1 ~ /^@/'
     ```
   - **Why MAPQ == 0 is Retained (Multi-mappers)**: Perfect multi-mappers (alignments mapping to multiple locations with equal score, mapped with MAPQ=0 by BWA) are *deliberately retained*. This is crucial because the IS element may insert into multicopy genes (such as ribosomal RNA operons, homologous transposon regions, or duplicated genomic segments). Discarding MAPQ=0 reads would completely hide insertions in these repetitive regions.
   - **Filtering Low-Quality Alignments (0 < MAPQ < min_mapq)**: Reads with intermediate or low quality mapping (where BWA favors one position slightly but with low confidence) are discarded to avoid false positive calls arising from partial/spurious sequence homologies.
   - **High-Confidence Alignments (MAPQ >= min_mapq)**: Standard unique alignments (default threshold MAPQ >= 30) are preserved.
3. **Peak Calling and Coverage Depth Cutoff**:
   - The depth of mapped reads is calculated genome-wide independently for left (5'/HEAD) and right (3'/TAIL) flanks.
   - To filter out low-coverage background noise, a minimal coverage filter is applied:
     ```bash
     awk '$4 >= {cutoff}' {left_cov} > {left_final_cov}
     ```
   - Only bases with a coverage depth greater than or equal to the `--cutoff` threshold (default 6) are preserved.
   - The remaining coordinates are merged using `bedtools merge -d {merging}` to call candidate flanking peaks. Any peak analyzed by the pairing logic is therefore supported by at least 6 reads.

### 4.2 Chromosome and Plasmid References
Bacteria often harbor plasmids in addition to the primary chromosome.
1. **Combined Reference Requirement**:
   - If an IS element is located on a plasmid, but only the chromosome is provided as the reference sequence, plasmid-derived reads may remain unmapped or mistakenly align with low-confidence to chromosomal regions sharing weak sequence homology.
   - To overcome this, a combined reference file (in a multi-fasta or multi-genbank file) containing both the chromosome and all plasmids can be used.
2. **Simultaneous Mapping (Parallel Competition)**:
   - When a combined reference containing both the chromosome and plasmids is provided, mapping is performed **simultaneously** rather than sequentially.
   - ISLAN parses all records from the input file, merges them into a single temporary reference FASTA, and builds a single combined BWA index.
   - During `bwa mem` execution, all reads are aligned against this combined index in a single run. This allows reads to compete across all replicons (chromosome and plasmids) at the same time. This is critical for evaluating correct mapping quality (MAPQ) and identifying true multi-mappers vs. unique plasmid/chromosomal insertions.
3. **Ambiguous Insertions (Similar/Duplicated Regions)**:
   - When identical sequence segments or homologous genes are present on both the chromosome and a plasmid, or duplicated on the chromosome:
     - The flanking reads will map to all copies with MAPQ=0.
     - The pipeline will report these as parallel hits on all homologous regions, reflecting the mathematical ambiguity of the insertion.
     - Users can resolve these by checking coverage depth (plasmids usually have higher copy numbers and thus higher coverage than the chromosome) or using long-read sequencing verification.

### 4.3 The Three-Stage Pairing Logic

#### Step 0: Chimera Filtering (Pre-filtering)
Before executing the three-stage pairing stages, the algorithm performs PCR chimera filtering on the raw peak pools:
- **Condition**: If a `HEAD` (left) and `TAIL` (right) peak overlap by 90%+ (`is_full_overlap`) and their coverage depth ratio exceeds `5.0`.
- **Action**: The minor peak is flagged as a chimera and immediately discarded from the active pools. It is excluded from all subsequent Stage 1, 2, or 3 pairing/singleton checks.

#### Stage 1: Resolve Known (Endogenous) IS Elements
1. **Target-Guided Scan**:
   - The reference genome is mapped against target elements carrying the `:FULL` suffix in the [targets.fasta](file:///data/SiNguyen/1.SIXTEEN/IS-SEQ/PISE/config/targets.fasta) file using `bwa mem -a` (or BLASTN).
   - This scan determines the exact genomic coordinates and orientation of all known target copies present in the reference sequence.
2. **Peak Pairing**:
   - For each known IS element coordinate interval `[start, end]` and orientation:
     - In **Forward (`+`)** orientation: The algorithm searches for a `HEAD` (5') peak near `start`, and a `TAIL` (3') peak near `end` within a window of `--flank-len`.
     - In **Reverse (`-`)** orientation: The algorithm searches for a `TAIL` (3') peak near `start`, and a `HEAD` (5') peak near `end`.
   - If both peaks are detected, they are paired as a `Known Pair`, marked as `paired`, and removed from the active pools.
3. **Logging**: The pipeline reports the total known copies found on the reference and how many were successfully paired.

#### Stage 2: Resolve Novel Insertions and Tandems
Active un-paired peaks are processed using coordinate-based overlap check functions:
1. **Novel Same-Direction Tandem**:
   - Signature: Two insertions in the same direction (`++` or `--`) next to each other.
   - Coordinate check: `l_peak_1` (`HEAD`) partially overlaps with a central pair consisting of `r_peak_1` (`TAIL`) and `l_peak_2` (`HEAD`) that fully overlap, which in turn partially overlaps with `r_peak_2` (`TAIL`).
   - Reported as `Tandem Gap (same direction)`.
2. **Novel Opposite-Direction Tandem**:
   - Signature `+-` (pointing towards each other): Two `HEAD` (5') peaks (`lp1`, `lp2`) both partially overlap with a single central `rpN` (`TAIL`). Reported as `Tandem Pair (+-)`.
   - Signature `-+` (pointing away from each other): Two `TAIL` (3') peaks (`rp1`, `rp2`) both partially overlap with a single central `lpN` (`HEAD`). Reported as `Tandem Pair (-+)`.
3. **Novel Pair (with Target Site Duplications - TSD)**:
   - Signature: A remaining `HEAD` (5') and `TAIL` (3') peak partially overlap each other.
   - Distance logic: The overlap size must be smaller than the `MAX_PAIRING_DISTANCE` (default: 100 bp).
   - Flagging possible false positives: If the overlap size is larger than the `MAX_TSD_OVERLAP` threshold (default: 20 bp), the call is appended with a `*` suffix (i.e. `novel (TSD)*`) to pinpoint potential empty/wild-type loci arising from non-specific primer binding. Note that this suffix is only biologically meaningful for and restricted to the `novel (TSD)` class.
   - Reported as `Novel Pair (TSD)` (or `Novel Pair (TSD)*`).
4. **Novel Pair (Standard)**:
   - Signature: A remaining `HEAD` (5') and `TAIL` (3') peak do not overlap but are within `MAX_PAIRING_DISTANCE` (default: 100 bp) of each other.
   - Reported as `Novel Pair`.

#### Stage 3: Resolve Singletons and Noise
1. **Full Flank Overlap**:
    - If a remaining `HEAD` and `TAIL` peak fully overlap (containment ratio $> 80\%$), they represent co-locating flank coverage.
    - **Biological & Technical Sources**:
      - **Transposition Intermediates (Hairpins / IS-Circles)**: Copy-out / paste-in transposition (e.g. in IS3, IS30, IS256 families) forms single-stranded hairpin or circularized IS intermediates joining IR-L and IR-R prior to target insertion, generating symmetric overlapping read flanks.
      - **Off-Target PCR Amplification & Chimeras**: Non-specific primer annealing or chimeric PCR extension fragments covering a non-target region.
    - Classified as `Full Flank Overlap` in result tables and report summaries (moved to `_unpaired.tsv`).
2. **Singletons**:
    - Remaining un-paired peaks are reported as `HEAD-only` (5' flank only) or `TAIL-only` (3' flank only) singletons in the unpaired TSV.

---

## 5. Full Flank Overlap Classification

Full Flank Overlap calls are identified purely based on the overlap fraction between co-locating peaks:
*   **Rule**: If a `HEAD` and `TAIL` peak on the same chromosome overlap by 90%+ (`is_full_overlap`), the pair is classified as **`Full Flank Overlap`** (and output to `_unpaired.tsv`).
*   **No Depth Ratio Cutoff**: All fully overlapping peak pairs ($\ge 90\%$ containment) are retained as `Full Flank Overlap` calls without filtering by depth asymmetry ratios.
        style K fill:#e74c3c,stroke:#c0392b,stroke-width:2px,color:#fff
        style G fill:#f39c12,stroke:#d35400,stroke-width:2px,color:#fff
    end
```

### Resolving IS Orientation
Based on outward-facing primers, the orientation of a paired insertion is deduced using their coordinates:
*   **Forward (`F` / `+`) Orientation**: The 5' end (`HEAD`) is on the genomic left, and the 3' end (`TAIL`) is on the genomic right (`HEAD_coordinate < TAIL_coordinate`).
*   **Reverse (`R` / `-`) Orientation**: The 3' end (`TAIL`) is on the genomic left, and the 5' end (`HEAD`) is on the genomic right (`TAIL_coordinate < HEAD_coordinate`).

---

## 6. Output Data Format & Interpretation

ISLAN outputs two main TSV files:
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
| `call` | Classification of the hit: `known`, `novel`, `novel (TSD)`. If a `novel (TSD)` hit has a flanking overlap exceeding `MAX_TSD_OVERLAP` (20 bp), it is reported as `novel (TSD)*` to indicate possible false positives (empty/wild-type loci). Other classes do not receive the `*` suffix. |
| `left_pos` | Chromosome coordinate range of the left flanking region (derived from 5'/HEAD or 3'/TAIL depending on orientation). |
| `right_pos` | Chromosome coordinate range of the right flanking region (derived from 3'/TAIL or 5'/HEAD depending on orientation). |
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

---

## 7. HTML Visualization Report

The pipeline generates an interactive HTML visualization report (`{sample}__{reference}_report.html`) containing three main sections:

### 7.1 Section 1: Summary Table
- Displays a tight, compact summary table detailing the counts for each of the core detection classes (including zero-count rows shown in muted grey for easy follow).

### 7.2 Section 2: Known IS Loci on Reference (BLASTN)
- Displays all reference IS positions scanned by BLASTN (within a window defined by `EXTENSION_PADDING = 2500` bp).
- **POSITIVE (Green Badge)**: Indicates known IS loci with paired read evidence (plots only show the paired flanks and their mapped reads).
- **NEGATIVE (Amber Badge)**: Indicates known IS loci without paired read evidence (plots show singleton/unpaired flanking reads, GC content, gene annotations, and the IS body coordinates).

### 7.3 Section 3: Novel IS Loci
- Displays interactive alignments (using the `generate_combined_alignment_plotly` viewer) for each novel insertion.
- Cards are color-coded based on the detection class:
  - **Indigo border**: Standard `novel` insertions
  - **Teal border**: `novel (TSD)` insertions
  - **Amber border**: `novel (TSD)*` potential false-positive insertions
- CDS/gene disruption badges and flanking gene descriptions are displayed directly in the header of each card.
