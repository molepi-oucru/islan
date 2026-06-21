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

def process_batch_worker(batch_data, primers, min_cov, min_identity, min_len, index_len, min_len_forward, expected_index_seq=None, i5_mismatch=2, len_actual_primer=None):
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.gap_score = -1.0

    passed_filtered = []
    passed_head = []
    passed_tail = []

    for record_tuple in batch_data:
        header, read_seq, spacer, qual = record_tuple
        
        # 1. Length filter
        if len(read_seq) < min_len_forward:
            continue
            
        # 2. Expected index filter (only if filtering raw read input directly)
        if expected_index_seq is not None:
            idx_seq = read_seq[:index_len]
            if sum(c1 != c2 for c1, c2 in zip(idx_seq, expected_index_seq)) > i5_mismatch:
                continue

        # 3. Trim i5 index
        trimmed_seq = read_seq[index_len:]
        trimmed_qual = qual[index_len:]
        trimmed_record_str = f"{header}{trimmed_seq}\n{spacer}{trimmed_qual}\n"
        
        passed_filtered.append(trimmed_record_str)

        # 4. Local alignment to HEAD/TAIL
        best_match = None
        for p_type in ['HEAD', 'TAIL']:
            if p_type not in primers: continue
            p_seq = primers[p_type]
            
            read_head = trimmed_seq[:len(p_seq)]
            alignments = aligner.align(p_seq, read_head)
            if alignments:
                alignment = alignments[0]
                target_start = alignment.aligned[0][0][0]
                target_end = alignment.aligned[0][-1][1]
                cov = (target_end - target_start) / len(p_seq)
                
                counts = alignment.counts()
                align_len = counts.identities + counts.mismatches + counts.gaps
                identity = counts.identities / align_len if align_len > 0 else 0.0
                
                is_std_match = (cov >= min_cov) and (identity >= min_identity)
                is_adv_match = False
                
                if len_actual_primer is not None and not is_std_match:
                    if len(trimmed_seq) >= len_actual_primer:
                        read_sub = trimmed_seq[:len_actual_primer]
                        primer_sub = p_seq[:len_actual_primer]
                        if read_sub == primer_sub:
                            is_adv_match = True
                            
                if is_std_match or is_adv_match:
                    score = alignment.score if is_std_match else float(len_actual_primer)
                    if best_match is None or score > best_match['score']:
                        best_match = {'type': p_type, 'score': score, 'primer_len': len(p_seq)}
                        
        if best_match is not None:
            p_len = best_match['primer_len']
            extracted_seq = trimmed_seq[:p_len]
            extracted_qual = trimmed_qual[:p_len]
            extracted_record_str = f"{header}{extracted_seq}\n{spacer}{extracted_qual}\n"
            if best_match['type'] == 'HEAD':
                passed_head.append(extracted_record_str)
            elif best_match['type'] == 'TAIL':
                passed_tail.append(extracted_record_str)
                
    return passed_filtered, passed_head, passed_tail, len(batch_data)

def stream_batches(fastq_path, batch_size=5000):
    with gzip.open(fastq_path, "rt") as infile:
        current_batch = []
        while True:
            header = infile.readline()
            if not header: break
            seq = infile.readline().strip()
            spacer = infile.readline()
            qual = infile.readline().strip()
            
            # Skip poly-N reads (already dropped by QC, but as fallback)
            if seq and seq[0] == 'N' and seq == 'N' * len(seq): continue
                
            current_batch.append((header, seq, spacer, qual))
            
            if len(current_batch) == batch_size:
                yield current_batch
                current_batch = []
        if current_batch:
            yield current_batch

def process_forward_reads(forward_fastq_path, primers, min_cov, min_identity, min_len, index_len, min_len_forward,
                          output_filtered_path, output_head_path, output_tail_path,
                          expected_index_seq=None, i5_mismatch=2, len_actual_primer=None, threads=4):
    logging.info(f"Filtering forward reads from: {forward_fastq_path}")
    logging.info(f"Targeting outputs: {output_filtered_path}, {output_head_path}, {output_tail_path}")
    
    os.makedirs(os.path.dirname(output_filtered_path), exist_ok=True)
    
    total_input = 0
    total_filtered = 0
    total_head = 0
    total_tail = 0

    worker_func = partial(
        process_batch_worker,
        primers=primers,
        min_cov=min_cov,
        min_identity=min_identity,
        min_len=min_len,
        index_len=index_len,
        min_len_forward=min_len_forward,
        expected_index_seq=expected_index_seq,
        i5_mismatch=i5_mismatch,
        len_actual_primer=len_actual_primer
    )

    try:
        with gzip.open(output_filtered_path, "wt") as out_filt, \
             gzip.open(output_head_path, "wt") as out_head, \
             gzip.open(output_tail_path, "wt") as out_tail:
             
            batches_gen = stream_batches(forward_fastq_path, batch_size=5000)

            with multiprocessing.Pool(processes=threads) as pool:
                for passed_filtered, passed_head, passed_tail, batch_size in pool.imap(worker_func, batches_gen):
                    for record_str in passed_filtered:
                        out_filt.write(record_str)
                    for record_str in passed_head:
                        out_head.write(record_str)
                    for record_str in passed_tail:
                        out_tail.write(record_str)

                    total_filtered += len(passed_filtered)
                    total_head += len(passed_head)
                    total_tail += len(passed_tail)
                    total_input += batch_size

    except Exception as e:
        logging.error(f"Error during parallel filtering: {e}")
        sys.exit(1)

    passed_pct = (total_filtered / total_input * 100) if total_input > 0 else 0
    head_pct = (total_head / total_filtered * 100) if total_filtered > 0 else 0
    tail_pct = (total_tail / total_filtered * 100) if total_filtered > 0 else 0

    logging.info(f"Filtering Complete: {total_filtered} reads passed length filter ({passed_pct:.2f}% of input reads).")
    logging.info(f"  - Classified HEAD: {total_head} ({head_pct:.2f}%)")
    logging.info(f"  - Classified TAIL: {total_tail} ({tail_pct:.2f}%)")

    return {
        'total_filtered': total_filtered,
        'total_head': total_head,
        'total_tail': total_tail
    }

def main():
    pass
