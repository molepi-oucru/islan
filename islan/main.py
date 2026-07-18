import os
import sys
import argparse
import logging
from islan.constants import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_MIN_CLIP,
    DEFAULT_MAX_CLIP,
    DEFAULT_CUTOFF,
    DEFAULT_MERGING,
    DEFAULT_IS_LENGTH,
    DEFAULT_MIN_MAPQ,
    DEFAULT_FLANK_LEN,
    DEFAULT_THREADS
)

def str2bool(v):
    if isinstance(v, bool): return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'): return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'): return False
    else: raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    parser = argparse.ArgumentParser(description="islan: Insertion Sequence Landscape Analyzer for IS-Seq sequencing data analysis.")
    subparsers = parser.add_subparsers(dest="subcommand", required=True, help="Subcommand to run")

    # Subparser: pre-process
    preprocess_parser = subparsers.add_parser("pre-process", aliases=["preprocess"], help="QC, demultiplexing, and filtering on forward reads")
    preprocess_parser.add_argument("forward_reads", help="Path to forward reads FASTQ file (_1.fastq.gz) or raw reads directory")
    preprocess_parser.add_argument("reverse_reads", nargs='?', default=None, help="Path to reverse reads FASTQ file (optional)")
    preprocess_parser.add_argument("-c", "--config", default=None, help="Path to pipeline configuration YAML file")
    preprocess_parser.add_argument("-t", "--threads", type=int, default=None, help="Number of threads/processes")
    preprocess_parser.add_argument("--index_i5", type=str2bool, default=None, help="Enable 8 bp i5 index prefix sequence offset")
    preprocess_parser.add_argument("--qc", type=str2bool, default=None, help="Run QC and demultiplexing step (default: True)")
    preprocess_parser.add_argument("--i5-mismatch", "--i5_mismatch", type=int, default=None, help="Mismatch tolerance for i5 index demultiplexing")
    preprocess_parser.add_argument("--min_len", type=int, default=None, help="Minimum remaining length threshold for primer matching")
    preprocess_parser.add_argument("-o", "--output_dir", default=None, help="Main pipeline output directory")

    # Subparser: asv-analysis
    asv_parser = subparsers.add_parser("asv-analysis", help="Run ASV analysis (placeholder)")
    asv_parser.add_argument("-c", "--config", default=None, help="Path to pipeline configuration YAML file")

    # Subparser: pairing
    pairing_parser = subparsers.add_parser("pairing", help="Extract reverse reads matching forward reads or ID list")
    pairing_parser.add_argument("-f", "--forward", required=True, help="Filtered forward FASTQ (.gz) or ID text file")
    pairing_parser.add_argument("-r", "--reverse", required=True, help="Raw reverse FASTQ (.gz)")
    pairing_parser.add_argument("-o", "--output", required=True, help="Output reverse FASTQ (.gz)")

    # Subparser: is-mapping
    ismapper_parser = subparsers.add_parser("is-mapping", help="Identify IS insertion sites (adapted from ISMapper)")
    ismapper_parser.add_argument("-c", "--config", default=None, help="Path to pipeline config YAML. CLI args override config values.")
    ismapper_parser.add_argument("--reads", nargs='+', help="Input reads (WGS mode)")
    ismapper_parser.add_argument("--queries", help="IS query FASTA (WGS mode)")
    ismapper_parser.add_argument("--reference", required=True, help="Reference genome FASTA/GenBank")
    ismapper_parser.add_argument("--targeted", action="store_true", help="Run in targeted mode (bypasses query mapping)")
    ismapper_parser.add_argument("--head", help="HEAD primer extracted reads (Targeted mode)")
    ismapper_parser.add_argument("--tail", help="TAIL primer extracted reads (Targeted mode)")
    ismapper_parser.add_argument("--filtered_forward", help="Filtered forward reads _1.fastq.gz (Targeted mode)")
    ismapper_parser.add_argument("--filtered_reverse", help="Filtered reverse reads _2.fastq.gz (Targeted mode)")
    ismapper_parser.add_argument("--forward-only", "--forward_only", dest="forward_only", action="store_true", help="Run analysis with forward reads only (single-end)")
    ismapper_parser.add_argument("-o", "--output_dir", default=DEFAULT_OUTPUT_DIR, help="Output directory")
    ismapper_parser.add_argument("-t", "--threads", type=int, default=None, help="Number of threads")
    ismapper_parser.add_argument("--min_clip", type=int, default=None, help="Minimum soft-clip size (WGS)")
    ismapper_parser.add_argument("--max_clip", type=int, default=None, help="Maximum soft-clip size (WGS)")
    ismapper_parser.add_argument("--cutoff", type=int, default=None, help="Minimum depth cutoff for reporting")
    ismapper_parser.add_argument("--merging", type=int, default=None, help="Merge distance for bedtools")
    ismapper_parser.add_argument("--is_length", type=int, default=None, help=f"Max gap distance to pair endogenous IS flanks (default {DEFAULT_IS_LENGTH})")
    ismapper_parser.add_argument("--min-mapq", "--min_mapq", dest="min_mapq", type=int, default=None, help="Minimum mapping quality to filter reads (retains MAPQ == 0 multi-mappers; discards 0 < MAPQ < min_mapq)")
    ismapper_parser.add_argument("--flank-len", "--flank_len", dest="flank_len", type=int, default=None, help=f"Flanking window around known IS element boundaries (default {DEFAULT_FLANK_LEN} bp)")
    ismapper_parser.add_argument("--temp", action="store_true", help="Keep the temporary files directory after successful completion (default: remove)")
    ismapper_parser.add_argument("--is_name", default=None, help="IS element name to filter known positions (e.g. ISKpn26). If not set, auto-derived from --queries filename (WGS mode) or all IS elements are used.")

    args = parser.parse_args()

    if args.subcommand == "asv-analysis":
        print("ASV-analysis subcommand placeholder. Detailed implementation not specified.")
        sys.exit(0)

    if args.subcommand == "pairing":
        from islan.extract_pairs import extract_pairs
        extract_pairs(args.forward, args.reverse, args.output)
        sys.exit(0)

    if args.subcommand in ("pre-process", "preprocess"):
        from islan.preprocess import run_preprocess
        run_preprocess(
            forward_reads=args.forward_reads,
            reverse_reads=args.reverse_reads,
            config_path=args.config,
            threads=args.threads,
            index_i5=args.index_i5,
            qc=args.qc,
            i5_mismatch=args.i5_mismatch,
            min_len=args.min_len,
            output_dir=args.output_dir
        )
        sys.exit(0)

    if args.subcommand == "is-mapping":
        from islan.mapping.main import run_is_mapping
        run_is_mapping(args)
        sys.exit(0)

if __name__ == "__main__":
    main()
