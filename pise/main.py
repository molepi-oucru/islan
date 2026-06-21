import os
import sys
import argparse

def str2bool(v):
    if isinstance(v, bool): return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'): return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'): return False
    else: raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(description="pise: Pipeline for IS-Seq sequencing data analysis.")
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

    args = parser.parse_args()

    if args.subcommand == "asv-analysis":
        print("ASV-analysis subcommand placeholder. Detailed implementation not specified.")
        sys.exit(0)

    if args.subcommand == "pairing":
        from pise.extract_pairs import extract_pairs
        extract_pairs(args.forward, args.reverse, args.output)
        sys.exit(0)

    if args.subcommand in ("pre-process", "preprocess"):
        from pise.preprocess import run_preprocess
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

if __name__ == "__main__":
    main()
