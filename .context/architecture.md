# Pipeline Architecture & Data Flow

This document details the modular architecture of the IS-Seq analysis pipeline.

## Modular Component Diagram

```mermaid
graph TD
    A1[Raw Forward FASTQ Reads _1] --> B(qc.py / Demux)
    C[config/primers.fasta] --> B
    B --> D[Demultiplexed Forward Reads]
    D --> B1(filter_reads.py)
    B1 --> D2[Filtered Forward Reads]
    
    A2[Raw Reverse FASTQ Reads _2] --> B2(extract_pairs.py)
    D2 --> B2
    B2 --> D3[Filtered Reverse Reads]
    
    D2 --> E(Alignment Module - bwa-mem2)
    D3 --> E
    F[Reference Genome] --> E
    E --> G[Mapped BAM File]
    G --> H(Junction Caller)
    H --> I[Quantified Insertion Sites & Reports]
```

## Configuration Schema (`config.yaml`)
The pipeline parameterizes its modules via `config/config.yaml`. The schema is defined as:

```yaml
# Global Parameters
target_is_element: "IS1R_IS1"
primers_file: "config/primers.fasta"
output_dir: "results"

# Preprocessing Step Parameters
preprocessing:
  min_cov: 0.98
  min_identity: 0.9
  min_len: 20
  index_i5: true
  len_actual_primer: null # Optional, set to integer to enable advanced rescue
  threads: 4 # Number of parallel processes to use for read filtering
  qc: true # Enable/disable QC and demultiplexing step (default: true)

# Alignment Step Parameters (Future)
alignment:
  reference_genome: "data/reference/genome.fasta"
  threads: 4
```
