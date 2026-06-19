import gzip
import os
import argparse
import logging
import sys

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def extract_pairs(forward_fastq, reverse_fastq_raw, output_reverse_fastq):
    """
    Extracts reads from reverse_fastq_raw that have an ID present in forward_fastq.
    This uses native 4-line parsing for extreme speed.
    """
    logging.info(f"Extracting valid pairs from {reverse_fastq_raw} based on {forward_fastq}")
    
    valid_ids = set()
    # 1. Collect valid IDs from forward file
    try:
        with gzip.open(forward_fastq, "rt") as f_in:
            while True:
                header = f_in.readline()
                if not header: break
                f_in.readline() # seq
                f_in.readline() # spacer
                f_in.readline() # qual
                # Extract ID: @M01234:56... 1:N:0:1 -> M01234:56...
                rec_id = header.strip()[1:].split()[0]
                valid_ids.add(rec_id)
    except FileNotFoundError:
        logging.error(f"Forward file not found: {forward_fastq}")
        sys.exit(1)
        
    logging.info(f"Collected {len(valid_ids)} valid read IDs from forward file.")
    
    # 2. Extract matching reverse reads
    extracted = 0
    total = 0
    try:
        with gzip.open(reverse_fastq_raw, "rt") as r_in, gzip.open(output_reverse_fastq, "wt") as r_out:
            while True:
                header = r_in.readline()
                if not header: break
                seq = r_in.readline()
                spacer = r_in.readline()
                qual = r_in.readline()
                
                total += 1
                rec_id = header.strip()[1:].split()[0]
                if rec_id in valid_ids:
                    r_out.write(f"{header}{seq}{spacer}{qual}")
                    extracted += 1
    except FileNotFoundError:
        logging.error(f"Reverse raw file not found: {reverse_fastq_raw}")
        sys.exit(1)
        
    logging.info(f"Extraction complete: {extracted}/{total} reverse reads kept.")
    logging.info(f"Output saved to: {output_reverse_fastq}")

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Extract reverse reads matching a filtered forward FASTQ")
    parser.add_argument("-f", "--forward", required=True, help="Filtered forward FASTQ (.gz)")
    parser.add_argument("-r", "--reverse", required=True, help="Raw reverse FASTQ (.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output reverse FASTQ (.gz)")
    args = parser.parse_args()
    
    extract_pairs(args.forward, args.reverse, args.output)

if __name__ == "__main__":
    main()
