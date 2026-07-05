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

def load_ids(file_path):
    import gzip
    ids = set()
    is_gz = file_path.endswith(".gz")
    open_func = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"
    
    try:
        with open_func(file_path, mode) as f:
            first_line = f.readline()
            if not first_line:
                return ids
                
            f.seek(0)
            if first_line.startswith("@"):
                # FASTQ format
                while True:
                    header = f.readline()
                    if not header: break
                    f.readline() # seq
                    f.readline() # spacer
                    f.readline() # qual
                    rec_id = header.strip()[1:].split()[0]
                    ids.add(rec_id)
            else:
                # Plain text format (one ID per line)
                for line in f:
                    stripped = line.strip()
                    if stripped:
                        ids.add(stripped.split()[0])
    except Exception as e:
        logging.error(f"Error reading ID file '{file_path}': {e}")
        sys.exit(1)
    return ids

def extract_pairs(id_or_fastq_path, reverse_fastq_raw, output_reverse_fastq):
    """
    Extracts reads from reverse_fastq_raw that have an ID present in id_or_fastq_path.
    Uses pyfastx for ultra-fast extraction.
    """
    import pyfastx
    logging.info(f"Extracting valid pairs from {reverse_fastq_raw} based on {id_or_fastq_path} using pyfastx")
    
    valid_ids = load_ids(id_or_fastq_path)
    logging.info(f"Loaded {len(valid_ids)} unique read IDs.")
    
    extracted = 0
    try:
        logging.info(f"Indexing reverse reads: {reverse_fastq_raw}")
        fq_rev = pyfastx.Fastx(reverse_fastq_raw)
        
        with gzip.open(output_reverse_fastq, "wt") as r_out:
            for name, seq, qual in fq_rev:
                name_base = name.split()[0]
                if name_base in valid_ids:
                    r_out.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    extracted += 1
                    
        logging.info(f"Extraction complete: {extracted} reverse reads kept.")
        logging.info(f"Output saved to: {output_reverse_fastq}")
    except Exception as e:
        logging.error(f"Error during extraction: {e}")
        sys.exit(1)

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Extract reverse reads matching a filtered forward FASTQ or ID list")
    parser.add_argument("-f", "--forward", required=True, help="Filtered forward FASTQ (.gz) or ID text file")
    parser.add_argument("-r", "--reverse", required=True, help="Raw reverse FASTQ (.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output reverse FASTQ (.gz)")
    args = parser.parse_args()
    
    extract_pairs(args.forward, args.reverse, args.output)

if __name__ == "__main__":
    main()
