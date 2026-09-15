# Adapted from ISMapper
# Copyright (c) 2014, Jane Hawkey, Kathryn Holt

import os
import re
import shlex
import logging
import subprocess
import gzip
import pyfastx
from Bio.Seq import Seq
from Bio import SeqIO
import pysam

def check_command(cmd_name):
    """Check if command exists and return the actual command to use (bwa-mem2 vs bwa)"""
    import shutil
    if shutil.which(cmd_name):
        return cmd_name
    return None

def run_command(cmd_list, shell=False):
    """Run a command using subprocess."""
    cmd_str = " ".join(cmd_list) if not shell else cmd_list
    logging.debug(f"Running: {cmd_str}")
    try:
        if shell:
            subprocess.run(cmd_list, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            subprocess.run(cmd_list, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd_str}")
        logging.error(f"Error output: {e.stderr.decode('utf-8')}")
        raise

def bwa_index(fasta):
    """Build bwa index for the fasta file using bwa-mem2 if available, else bwa."""
    bwa_cmd = check_command('bwa-mem2') or check_command('bwa')
    if not bwa_cmd:
        raise EnvironmentError("Neither bwa-mem2 nor bwa is installed.")
    
    # check if index exists
    if bwa_cmd == 'bwa-mem2':
        built_index = fasta + '.bwt.2bit.64'
    else:
        built_index = fasta + '.bwt'
        
    if os.path.exists(built_index):
        logging.info(f'Index for {fasta} already exists.')
    else:
        logging.info(f'Building {bwa_cmd} index for {fasta}...')
        run_command([bwa_cmd, 'index', fasta])
    return bwa_cmd



def get_ids_and_primer_len(fastq_file):
    """Extract IDs and the primer length from a HEAD or TAIL file."""
    ids = []
    primer_len = 0
    if not os.path.exists(fastq_file):
        return ids, primer_len
        
    try:
        with gzip.open(fastq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    rec_id = line.strip()[1:].split()[0]
                    ids.append(rec_id)
                elif i % 4 == 1 and primer_len == 0:
                    primer_len = len(line.strip())
    except EOFError:
        pass
    return ids, primer_len

def extract_targeted_flanks_pyfastx(filtered_1_fastq, filtered_2_fastq, head_fastq, tail_fastq, tmp_folder, sample_prefix, forward_only=False):
    """
    Targeted Mode: Use pyfastx to rapidly extract flanks from targeted IS-Seq data.
    """
    left_final = os.path.join(tmp_folder, f"{sample_prefix}_left_final.fastq")
    right_final = os.path.join(tmp_folder, f"{sample_prefix}_right_final.fastq")
    
    head_ids, head_primer_len = get_ids_and_primer_len(head_fastq)
    tail_ids, tail_primer_len = get_ids_and_primer_len(tail_fastq)
    
    if not head_ids:
        logging.warning(f"HEAD file {head_fastq} is empty or missing.")
    if not tail_ids:
        logging.warning(f"TAIL file {tail_fastq} is empty or missing.")
        
    head_id_set = set(head_ids)
    tail_id_set = set(tail_ids)
    
    logging.info(f"Loaded {len(head_id_set)} HEAD IDs and {len(tail_id_set)} TAIL IDs.")
    
    # 1. Build index and extract forward reads (flanks)
    logging.info(f"Indexing forward reads: {filtered_1_fastq}")
    fq_fwd = pyfastx.Fastx(filtered_1_fastq)
    
    with open(left_final, 'w') as out_left, open(right_final, 'w') as out_right:
        for name, seq, qual in fq_fwd:
            name_base = name.split()[0]
            if name_base in head_id_set:
                flank_seq = seq[head_primer_len:]
                flank_qual = qual[head_primer_len:]
                if len(flank_seq) > 0:
                    out_left.write(f"@{name}\n{flank_seq}\n+\n{flank_qual}\n")
            elif name_base in tail_id_set:
                flank_seq = seq[tail_primer_len:]
                flank_qual = qual[tail_primer_len:]
                if len(flank_seq) > 0:
                    out_right.write(f"@{name}\n{flank_seq}\n+\n{flank_qual}\n")
                    
    if not forward_only:
        # 2. Extract matching reverse reads
        logging.info(f"Indexing reverse reads: {filtered_2_fastq}")
        fq_rev = pyfastx.Fastx(filtered_2_fastq)
        
        left_rev = os.path.join(tmp_folder, f"{sample_prefix}_left_rev.fastq")
        right_rev = os.path.join(tmp_folder, f"{sample_prefix}_right_rev.fastq")
        
        with open(left_rev, 'w') as out_left, open(right_rev, 'w') as out_right:
            for name, seq, qual in fq_rev:
                name_base = name.split()[0]
                if name_base in head_id_set:
                    out_left.write(f"@{name}\n{seq}\n+\n{qual}\n")
                elif name_base in tail_id_set:
                    out_right.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    
        # Interleave them for BWA PE mapping if necessary, or just treat as single end. 
        # ISMapper treats extracted unmapped reads as single end for mapping to ref.
        # Wait, ISMapper treats the flanking reads as single-end in map_to_ref_seq:
        # run_command(['bwa', 'mem', '-t', bwa_threads, ref_seq_file, left_flanking, ...])
        # So we don't necessarily need the reverse reads here unless we map paired.
        # We will cat them together to mimic ISMapper's behavior, which pools all left reads.
        
        # Mimic ISMapper pooling:
        run_command(f"cat {left_rev} >> {left_final}", shell=True)
        run_command(f"cat {right_rev} >> {right_final}", shell=True)
    
    return left_final, right_final

def map_to_ref_seq(ref_fasta, sample_prefix, left_flanking, right_flanking, tmp_folder, out_folder, threads, min_mapq=30):
    """
    Map the extracted HEAD and TAIL reads to the reference genome independently.
    """
    bwa_cmd = bwa_index(ref_fasta)
    ref_base = os.path.basename(ref_fasta).rsplit('.', 1)[0]
    
    left_sorted = os.path.join(out_folder, f"{sample_prefix}_left_{ref_base}.sorted.bam")
    right_sorted = os.path.join(out_folder, f"{sample_prefix}_right_{ref_base}.sorted.bam")
    
    import subprocess
    is_legacy_samtools = False
    try:
        res = subprocess.run(["samtools", "sort"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = res.stderr or res.stdout
        if "-o FILE" not in output and "output file name" not in output:
            is_legacy_samtools = True
    except Exception:
        pass

    # Ensure reference FASTA is indexed with faidx for BAM header safety
    if not os.path.exists(f"{ref_fasta}.fai"):
        try:
            run_command(f"samtools faidx {ref_fasta}", shell=True)
        except Exception:
            pass
            
    fai_opt = f"-t {ref_fasta}.fai" if os.path.exists(f"{ref_fasta}.fai") else ""

    if is_legacy_samtools:
        left_unsorted = left_sorted.replace(".sorted.bam", "_unsorted.bam")
        right_unsorted = right_sorted.replace(".sorted.bam", "_unsorted.bam")
        left_prefix = left_sorted.replace(".sorted.bam", "")
        right_prefix = right_sorted.replace(".sorted.bam", "")
        left_cmd = f"{bwa_cmd} mem -t {threads} {ref_fasta} {left_flanking} | awk '$5 == 0 || $5 >= {min_mapq} || $1 ~ /^@/' | samtools view -bS - > {left_unsorted} && samtools sort -@ {threads} {left_unsorted} {left_prefix} && mv {left_prefix}.bam {left_sorted} && rm -f {left_unsorted}"
        right_cmd = f"{bwa_cmd} mem -t {threads} {ref_fasta} {right_flanking} | awk '$5 == 0 || $5 >= {min_mapq} || $1 ~ /^@/' | samtools view -bS - > {right_unsorted} && samtools sort -@ {threads} {right_unsorted} {right_prefix} && mv {right_prefix}.bam {right_sorted} && rm -f {right_unsorted}"
    else:
        left_cmd = f"{bwa_cmd} mem -t {threads} {ref_fasta} {left_flanking} | samtools view -h | awk '$5 == 0 || $5 >= {min_mapq} || $1 ~ /^@/' | samtools sort -@ {threads} -o {left_sorted} -"
        right_cmd = f"{bwa_cmd} mem -t {threads} {ref_fasta} {right_flanking} | samtools view -h | awk '$5 == 0 || $5 >= {min_mapq} || $1 ~ /^@/' | samtools sort -@ {threads} -o {right_sorted} -"

    def _run_or_create_empty_bam(fastq_path, sorted_bam_path, ref_fasta_path, cmd_str):
        if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
            logging.info(f"Flanking FASTQ {fastq_path} is empty. Creating empty BAM file with reference header...")
            try:
                ref_records = list(SeqIO.parse(ref_fasta_path, "fasta"))
                header = {
                    'HD': {'VN': '1.0', 'SO': 'coordinate'},
                    'SQ': [{'SN': rec.id, 'LN': len(rec.seq)} for rec in ref_records]
                }
                with pysam.AlignmentFile(sorted_bam_path, "wb", header=header) as empty_bam:
                    pass
            except Exception as e:
                logging.warning(f"Could not write empty BAM with pysam: {e}. Running fallback command.")
                run_command(cmd_str, shell=True)
        else:
            run_command(cmd_str, shell=True)

    logging.info(f"Mapping left flanks to reference and sorting...")
    _run_or_create_empty_bam(left_flanking, left_sorted, ref_fasta, left_cmd)
    logging.info(f"Mapping right flanks to reference and sorting...")
    _run_or_create_empty_bam(right_flanking, right_sorted, ref_fasta, right_cmd)
    
    # Index BAMs
    logging.info(f"Indexing BAM files...")
    if is_legacy_samtools:
        run_command(f"samtools index {left_sorted}", shell=True)
        run_command(f"samtools index {right_sorted}", shell=True)
    else:
        run_command(f"samtools index -@ {threads} {left_sorted}", shell=True)
        run_command(f"samtools index -@ {threads} {right_sorted}", shell=True)
    
    return left_sorted, right_sorted

def create_bed_files(left_bam, right_bam, tmp_folder, out_folder, min_depth, merging):
    """
    Create Bedtools coverage maps and filter by depth cutoff.
    """
    sample_ref = os.path.basename(left_bam).replace(".sorted.bam", "")
    sample_prefix = sample_ref.split("_left_")[0]
    ref_base = sample_ref.split("_left_")[1]
    
    left_cov = os.path.join(tmp_folder, f"{sample_prefix}_left_{ref_base}_cov.bed")
    right_cov = os.path.join(tmp_folder, f"{sample_prefix}_right_{ref_base}_cov.bed")
    
    left_final_cov = os.path.join(out_folder, f"{sample_prefix}_left_{ref_base}_finalcov.bed")
    right_final_cov = os.path.join(out_folder, f"{sample_prefix}_right_{ref_base}_finalcov.bed")
    
    left_merged_bed = os.path.join(out_folder, f"{sample_prefix}_left_{ref_base}_merged.sorted.bed")
    right_merged_bed = os.path.join(out_folder, f"{sample_prefix}_right_{ref_base}_merged.sorted.bed")
    
    run_command(f"bedtools genomecov -ibam {left_bam} -bg > {left_cov}", shell=True)
    run_command(f"bedtools genomecov -ibam {right_bam} -bg > {right_cov}", shell=True)
    
    # Filter by min_depth using awk
    run_command(f"awk '$4 >= {min_depth}' {left_cov} > {left_final_cov}", shell=True)
    run_command(f"awk '$4 >= {min_depth}' {right_cov} > {right_final_cov}", shell=True)
    
    # Merge
    run_command(f"bedtools merge -d {merging} -i {left_final_cov} > {left_merged_bed}", shell=True)
    run_command(f"bedtools merge -d {merging} -i {right_final_cov} > {right_merged_bed}", shell=True)
        
    return left_merged_bed, right_merged_bed
