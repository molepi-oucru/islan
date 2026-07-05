#!/usr/bin/env python3
import os
import gzip
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Split a gzipped FASTQ file into HEAD and TAIL single-line FASTA files based on a sequence prefix pattern."
    )
    parser.add_argument(
        "fastq_gz", 
        help="Path to the input .fastq.gz file"
    )
    parser.add_argument(
        "pattern", 
        help="Prefix pattern for TAIL.fasta sequences (e.g., CTCTCTATTCAAAAT)"
    )
    parser.add_argument(
        "-o", "--output-dir", 
        default=".",
        help="Output directory for FASTA files (default: current directory)"
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.fastq_gz):
        print(f"Error: Input file '{args.fastq_gz}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
        
    # Determine basename and outputs
    basename = os.path.basename(args.fastq_gz)
    if basename.endswith(".fastq.gz"):
        basename_no_ext = basename[:-9]
    elif basename.endswith(".fq.gz"):
        basename_no_ext = basename[:-6]
    else:
        basename_no_ext = os.path.splitext(basename)[0]
        
    tail_path = os.path.join(args.output_dir, f"{basename_no_ext}.TAIL.fasta")
    head_path = os.path.join(args.output_dir, f"{basename_no_ext}.HEAD.fasta")
    
    print(f"Processing '{args.fastq_gz}'...")
    print(f"Pattern for TAIL: '{args.pattern}'")
    print(f"Output TAIL: {tail_path}")
    print(f"Output HEAD:   {head_path}")
    
    count_tail = 0
    count_head = 0
    
    try:
        with gzip.open(args.fastq_gz, "rt") as in_f:
            with open(tail_path, "w") as f_tail, open(head_path, "w") as f_head:
                while True:
                    header = in_f.readline()
                    if not header:
                        break
                    seq = in_f.readline().strip()
                    in_f.readline()  # + spacer
                    in_f.readline()  # quality
                    
                    # Parse sequence ID
                    seq_id = header.strip().split()[0][1:]
                    
                    fasta_record = f">{seq_id}\n{seq}\n"
                    if seq.startswith(args.pattern):
                        f_tail.write(fasta_record)
                        count_tail += 1
                    else:
                        f_head.write(fasta_record)
                        count_head += 1
                        
        print(f"Success! Processed {count_tail + count_head} reads:")
        print(f"  - TAIL reads: {count_tail}")
        print(f"  - HEAD reads:   {count_head}")
        
    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
