import os
import sys
import gzip
import logging
import multiprocessing
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq

def load_primers(fasta_file_path, target_is_element):
    primers = {}
    try:
        with open(fasta_file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split(':')
                if len(parts) >= 2 and parts[1] == target_is_element:
                    if ":HEAD" in record.id:
                        primers['HEAD'] = str(Seq(record.seq).reverse_complement())
                    elif ":TAIL" in record.id:
                        primers['TAIL'] = str(record.seq)
    except FileNotFoundError:
        logging.error(f"Primer FASTA file not found at {fasta_file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred while loading primers: {e}")
        sys.exit(1)
    return primers

def write_summary_tsv(filepath, stats):
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w') as f:
            f.write("Metric\tValue\n")
            for k, v in stats.items():
                f.write(f"{k}\t{v}\n")
        logging.info(f"Summary statistics saved to TSV: {filepath}")
    except Exception as e:
        logging.error(f"Failed to write summary statistics TSV: {e}")

def process_batch_worker(batch_data, primers, min_cov, min_identity, min_len, index_i5, len_actual_primer):
    index_len = 8 if index_i5 else 0
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.gap_score = -1.0

    passed_full = []
    passed_partial = []
    binned = []

    stats = {
        'passed_criterion1_full': 0,
        'passed_criterion1_partial': 0,
        'passed_criterion2_full': 0,
        'passed_criterion2_partial': 0,
        'head_matched_c1_full': 0,
        'tail_matched_c1_full': 0,
        'head_matched_c1_partial': 0,
        'tail_matched_c1_partial': 0,
        'head_matched_c2_full': 0,
        'tail_matched_c2_full': 0,
        'head_matched_c2_partial': 0,
        'tail_matched_c2_partial': 0
    }

    for record_tuple in batch_data:
        header, read_seq, spacer, qual = record_tuple
        
        # Poly-N reads are already dropped by QC, but we check just in case
        if all(c == 'N' for c in read_seq):
            continue

        best_match = None
        best_match_partial = None

        for p_type in ['HEAD', 'TAIL']:
            if p_type not in primers: continue
            p_seq = primers[p_type]

            if len(read_seq) >= index_len:
                read_head = read_seq[index_len : index_len + len(p_seq) + 5]
                alignments = aligner.align(p_seq, read_head)
                if alignments:
                    alignment = alignments[0]
                    target_start = alignment.aligned[0][0][0]
                    target_end = alignment.aligned[0][-1][1]
                    cov = (target_end - target_start) / len(p_seq)
                    
                    counts = alignment.counts()
                    align_len = counts.identities + counts.mismatches + counts.gaps
                    identity = counts.identities / align_len if align_len > 0 else 0.0

                    is_full_match = (cov >= min_cov) and (identity >= min_identity)
                    
                    if is_full_match:
                        score = alignment.score
                        if best_match is None or score > best_match['score']:
                            best_match = {
                                'type': p_type, 'score': score, 'matched_len': target_end,
                                'coverage': cov, 'identity': identity, 'is_full': True
                            }
                    elif len_actual_primer is not None:
                        if len(read_seq) - (index_len + target_end) >= min_len:
                            if len(read_seq) >= index_len + len_actual_primer:
                                read_sub = read_seq[index_len : index_len + len_actual_primer]
                                primer_sub = p_seq[:len_actual_primer]
                                if read_sub == primer_sub:
                                    score = float(len_actual_primer)
                                    if best_match_partial is None or score > best_match_partial['score']:
                                        best_match_partial = {
                                            'type': p_type, 'score': score, 'matched_len': target_end,
                                            'coverage': cov, 'identity': identity, 'is_full': False
                                        }

        matched_record = None
        if best_match is not None:
            matched_record = best_match
        elif best_match_partial is not None:
            matched_record = best_match_partial

        is_written = False
        out_str = f"{header}{read_seq}\n{spacer}{qual}\n"

        if matched_record is not None:
            is_full = matched_record['is_full']
            if is_full:
                stats['passed_criterion1_full'] += 1
                if matched_record['type'] == 'HEAD': stats['head_matched_c1_full'] += 1
                else: stats['tail_matched_c1_full'] += 1
            else:
                stats['passed_criterion1_partial'] += 1
                if matched_record['type'] == 'HEAD': stats['head_matched_c1_partial'] += 1
                else: stats['tail_matched_c1_partial'] += 1

            remaining_length = len(read_seq) - (index_len + matched_record['matched_len'])
            if remaining_length >= min_len:
                if is_full:
                    stats['passed_criterion2_full'] += 1
                    if matched_record['type'] == 'HEAD': stats['head_matched_c2_full'] += 1
                    else: stats['tail_matched_c2_full'] += 1
                    passed_full.append(out_str)
                    is_written = True
                else:
                    stats['passed_criterion2_partial'] += 1
                    if matched_record['type'] == 'HEAD': stats['head_matched_c2_partial'] += 1
                    else: stats['tail_matched_c2_partial'] += 1
                    passed_partial.append(out_str)
                    is_written = True

        if not is_written:
            binned.append(out_str)

    return passed_full, passed_partial, binned, stats

def stream_batches(fastq_path, batch_size=5000):
    with gzip.open(fastq_path, "rt") as infile:
        current_batch = []
        while True:
            header = infile.readline()
            if not header: break
            seq = infile.readline().strip()
            spacer = infile.readline()
            qual = infile.readline().strip()
            
            # Poly-N safety skip
            if all(c == 'N' for c in seq): continue
                
            current_batch.append((header, seq, spacer, qual))
            
            if len(current_batch) == batch_size:
                yield current_batch
                current_batch = []
        if current_batch:
            yield current_batch

def process_forward_reads(forward_fastq_path, primers, min_cov, min_identity, min_len, index_i5, len_actual_primer, output_forward_fastq_path, output_forward_partial_path=None, output_forward_bin_path=None, threads=4):
    logging.info(f"Processing forward reads from: {forward_fastq_path}")
    logging.info(f"Using {threads} threads/processes for filtering.")
    
    os.makedirs(os.path.dirname(output_forward_fastq_path), exist_ok=True)
    
    total_reads = 0
    passed_c1_full = passed_c1_partial = 0
    passed_c2_full = passed_c2_partial = 0
    head_matched_c1_full = tail_matched_c1_full = 0
    head_matched_c1_partial = tail_matched_c1_partial = 0
    head_matched_c2_full = tail_matched_c2_full = 0
    head_matched_c2_partial = tail_matched_c2_partial = 0

    worker_func = partial(
        process_batch_worker,
        primers=primers,
        min_cov=min_cov,
        min_identity=min_identity,
        min_len=min_len,
        index_i5=index_i5,
        len_actual_primer=len_actual_primer
    )

    try:
        with gzip.open(output_forward_fastq_path, "wt") as outfile_full, \
             gzip.open(output_forward_bin_path, "wt") as outfile_bin:
             
            outfile_partial = None
            if output_forward_partial_path:
                outfile_partial = gzip.open(output_forward_partial_path, "wt")

            batches_gen = stream_batches(forward_fastq_path, batch_size=5000)

            with multiprocessing.Pool(processes=threads) as pool:
                for passed_full, passed_partial, binned, batch_stats in pool.imap(worker_func, batches_gen):
                    for fq_str in passed_full:
                        outfile_full.write(fq_str)
                    if outfile_partial:
                        for fq_str in passed_partial:
                            outfile_partial.write(fq_str)
                    for fq_str in binned:
                        outfile_bin.write(fq_str)

                    # Update stats
                    total_reads += len(passed_full) + len(passed_partial) + len(binned)
                    passed_c1_full += batch_stats['passed_criterion1_full']
                    passed_c1_partial += batch_stats['passed_criterion1_partial']
                    passed_c2_full += batch_stats['passed_criterion2_full']
                    passed_c2_partial += batch_stats['passed_criterion2_partial']
                    
                    head_matched_c1_full += batch_stats['head_matched_c1_full']
                    tail_matched_c1_full += batch_stats['tail_matched_c1_full']
                    head_matched_c1_partial += batch_stats['head_matched_c1_partial']
                    tail_matched_c1_partial += batch_stats['tail_matched_c1_partial']
                    
                    head_matched_c2_full += batch_stats['head_matched_c2_full']
                    tail_matched_c2_full += batch_stats['tail_matched_c2_full']
                    head_matched_c2_partial += batch_stats['head_matched_c2_partial']
                    tail_matched_c2_partial += batch_stats['tail_matched_c2_partial']

            if outfile_partial:
                outfile_partial.close()

    except Exception as e:
        logging.error(f"Error during parallel filtering: {e}")
        sys.exit(1)

    return {
        'total_reads': total_reads,
        'passed_c1_full': passed_c1_full,
        'passed_c1_partial': passed_c1_partial,
        'passed_c2_full': passed_c2_full,
        'passed_c2_partial': passed_c2_partial,
        'head_c1_full': head_matched_c1_full,
        'tail_c1_full': tail_matched_c1_full,
        'head_c1_partial': head_matched_c1_partial,
        'tail_c1_partial': tail_matched_c1_partial,
        'head_c2_full': head_matched_c2_full,
        'tail_c2_full': tail_matched_c2_full,
        'head_c2_partial': head_matched_c2_partial,
        'tail_c2_partial': tail_matched_c2_partial
    }

def main():
    pass
