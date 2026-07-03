# Understanding Single-End vs. Paired-End Analysis in IS-Targeted Sequencing

This document explains the differences in insertion site detection (count, coordinates, and types) between forward-only (single-end) and paired-end analyses, and details the technical reasons behind these discrepancies.

## 1. Observed Differences

When comparing the analysis of forward-only reads versus paired-end reads in IS-targeted sequencing (like the `502-ISKpn26` dataset), three main differences are observed:

*   **Loss of Insertions & Increased "Unpaired" Hits**: Forward-only analysis often loses specific insertion calls (e.g., `novel (TSD)`). Instead of pairing the left and right flanks into a single insertion event, the algorithm categorizes them as separate, unlinked `HEAD` or `TAIL` fragments in the `unpaired.tsv` file.
*   **Coordinate Shifts**: The mapped coordinates of the flanks (the peaks) shift away from the true IS-genome junction. The gap between a `HEAD` and `TAIL` flank becomes artificially wider.
*   **Altered Insertion Types**: Because the gap between flanks widens, insertions that truly have a Target Site Duplication (TSD)—which typically appear as an overlap or a very short gap—are no longer recognized as TSDs. They are instead classified as standard `novel` insertions or dropped to `unpaired` if the gap exceeds the pairing threshold (e.g., > 100 bp).

## 2. Why do mapping depth and quality drop at the junction?

In targeted sequencing, the forward primer is designed on the IS element, pointing outward into the flanking genomic DNA. Therefore, **Read 1 (the forward read) starts exactly at, or very close to, the IS-genome junction.**

Despite being closest to the primer, mapping algorithms like BWA-MEM struggle precisely *because* the read starts at this complex boundary:

1.  **Seeding and Soft-Clipping**: BWA-MEM relies on finding Maximal Exact Matches (MEMs) to "seed" the alignment. The exact junction often contains non-genomic sequence (leftover primer, partial IS sequence) or biological complexities like microhomologies (TSDs). If the start of the read has mismatches or complex sequences, BWA-MEM cannot generate a high-quality seed there.
2.  **Algorithm Heuristics**: To prevent poor alignments, BWA-MEM applies a penalty for mismatches. If the start of the read drops the overall alignment score, the algorithm will perform **soft-clipping**—ignoring the first several bases of the read to maximize the score of the remaining genomic portion. 
3.  **Coordinate Shift**: Because the read is soft-clipped at the 5' end, the "start" of the mapped peak shifts downstream into the genome, away from the true junction. This artificially widens the measured distance between the left and right flanks of an insertion.

## 3. How do reverse reads help if they don't even reach the junction?

In a standard library preparation (e.g., Nextera tagmentation or mechanical shearing), the reverse read (Read 2) originates from a random fragmentation site downstream and reads *inward* toward the IS element. 

As you correctly noted, with a read length of ~150 bp and a fragment size of 300-500 bp, **Read 2 rarely reaches the actual IS junction.** However, it is still crucial for accurate mapping due to two major effects:

### A. Anchoring and Rescuing Read 1 (Mapping Quality)
Read 2 consists entirely of standard genomic DNA, far away from the messy IS junction. Consequently, BWA-MEM can map Read 2 with extremely high confidence (high MAPQ). 

In paired-end mode, BWA-MEM uses this high-quality Read 2 to "anchor" Read 1. The aligner knows that Read 1 must map near Read 2 and in the opposite orientation. This localized search constraint helps BWA-MEM overcome the seeding difficulties of Read 1 at the junction, forcing a more accurate alignment of Read 1 and reducing erroneous soft-clipping.

### B. Extending the Coverage Profile (Depth)
In forward-only mode, the mapping depth peak is strictly clamped to the length of Read 1 (~150 bp). If BWA-MEM soft-clips the start of the read, the entire 150 bp block shifts downstream.

In paired-end mode, the analysis tools use the **entire template fragment length** inferred from the read pair. The mapped region spans from the start of Read 1 to the end of Read 2 (e.g., 400 bp total). This creates a broad, highly stable coverage peak that physically connects the junction to the downstream genomic sequence. 
*   **Result**: The peak boundary remains firmly anchored at the true junction, the gap distance is accurately calculated (allowing detection of TSDs), and the left/right flanks are easily paired within the acceptable distance thresholds.

---

> [!NOTE] 
> **Summary**
> Single-end analysis relies entirely on the forward read, which struggles to map perfectly at the complex IS boundary, leading to soft-clipping, shifted coordinates, and widened gaps that break insertion pairing. Paired-end analysis uses the high-quality, downstream reverse read to "anchor" the alignment and extend the coverage profile, ensuring the true junction coordinates are preserved.
