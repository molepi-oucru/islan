import os
import sys
import logging
import argparse
from .mapping import map_to_is_query_wgs, extract_targeted_flanks_pyfastx, map_to_ref_seq, create_bed_files
from .reporting import create_typing_output
from .report_generator import generate_report
from islan.constants import DEFAULT_OUTPUT_DIR, MAPPING_LOG_FILE

def run_is_mapping(args):
    """Entry point for the is-mapping pipeline."""
    # Ensure directories
    out_dir = os.path.abspath(args.output_dir if args.output_dir else DEFAULT_OUTPUT_DIR)
    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Configure file logging
    log_file = os.path.join(out_dir, MAPPING_LOG_FILE)
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        if isinstance(handler, logging.FileHandler):
            root_logger.removeHandler(handler)
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    root_logger.addHandler(file_handler)
    
    sample_prefix = "sample" # Determine from input
    
    left_flanking = None
    right_flanking = None
    
    if args.targeted:
        logging.info("Running in TARGETED mode...")
        if args.forward_only:
            if not args.head or not args.tail or not args.filtered_forward:
                logging.error("Targeted forward-only mode requires --head, --tail, and --filtered_forward")
                sys.exit(1)
        else:
            if not args.head or not args.tail or not args.filtered_forward or not args.filtered_reverse:
                logging.error("Targeted mode requires --head, --tail, --filtered_forward, and --filtered_reverse")
                sys.exit(1)
            
        sample_prefix = os.path.basename(args.head).split("_")[0] if args.head else "sample"
        left_flanking, right_flanking = extract_targeted_flanks_pyfastx(
            filtered_1_fastq=args.filtered_forward,
            filtered_2_fastq=args.filtered_reverse,
            head_fastq=args.head,
            tail_fastq=args.tail,
            tmp_folder=tmp_dir,
            sample_prefix=sample_prefix,
            forward_only=args.forward_only
        )
    else:
        logging.info("Running in WGS mode...")
        # Assume reads are provided as two paired files for simplicity, or we can use the first two elements.
        forward_read = args.reads[0]
        if args.forward_only or len(args.reads) == 1:
            reverse_read = None
        else:
            reverse_read = args.reads[1]
        sample_prefix = os.path.basename(forward_read).split("_")[0]
        left_flanking, right_flanking = map_to_is_query_wgs(
            sample_prefix=sample_prefix,
            forward_fastq=forward_read,
            reverse_fastq=reverse_read,
            is_query_fasta=args.queries,
            tmp_folder=tmp_dir,
            out_dir=out_dir,
            min_clip=args.min_clip,
            max_clip=args.max_clip,
            threads=args.threads
        )
        
    logging.info("Flanking reads extracted. Proceeding to reference mapping.")
    
    # BWA requires FASTA format. If user provided a GenBank file, convert it to FASTA.
    ref_fasta = args.reference
    if ref_fasta.endswith('.gbk') or ref_fasta.endswith('.gb'):
        import Bio.SeqIO
        ref_base_name = os.path.basename(args.reference).rsplit('.', 1)[0]
        ref_fasta = os.path.join(tmp_dir, ref_base_name + '.fasta')
        if not os.path.exists(ref_fasta):
            logging.info(f"Converting GenBank to FASTA: {ref_fasta}")
            Bio.SeqIO.convert(args.reference, "genbank", ref_fasta, "fasta")
            
    # Locate targets or primers file
    targets_fasta = os.path.abspath("config/targets.fasta")
    if not os.path.exists(targets_fasta):
        targets_fasta = os.path.abspath("config/primers.fasta")

    left_bam, right_bam = map_to_ref_seq(
        ref_fasta=ref_fasta,
        sample_prefix=sample_prefix,
        left_flanking=left_flanking,
        right_flanking=right_flanking,
        tmp_folder=tmp_dir,
        out_folder=out_dir,
        threads=args.threads,
        min_mapq=args.min_mapq
    )
    
    left_merged, right_merged = create_bed_files(
        left_bam=left_bam,
        right_bam=right_bam,
        tmp_folder=tmp_dir,
        out_folder=out_dir,
        cutoff=args.cutoff,
        merging=args.merging
    )
    
    ref_base = os.path.basename(args.reference).rsplit('.', 1)[0]
    out_table = os.path.join(out_dir, f"{sample_prefix}__{ref_base}_table.tsv")
    
    create_typing_output(
        left_merged=left_merged,
        right_merged=right_merged,
        left_cov=os.path.join(tmp_dir, f"{sample_prefix}_left_{ref_base}_cov.bed"),
        right_cov=os.path.join(tmp_dir, f"{sample_prefix}_right_{ref_base}_cov.bed"),
        ref_fasta=args.reference,
        out_file=out_table,
        left_bam=left_bam,
        right_bam=right_bam,
        is_length=args.is_length,
        flank_len=args.flank_len,
        targets_fasta=targets_fasta,
        threads=args.threads
    )
    
    # Generate HTML report
    report_file = os.path.join(out_dir, f"{sample_prefix}__{ref_base}_report.html")
    generate_report(out_table, report_file, reference_file=args.reference, cutoff=args.cutoff)
    
    if not getattr(args, 'temp', False):
        import shutil
        logging.info("Cleaning up temporary files...")
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        logging.info(f"Keeping temporary files in {tmp_dir}")
        
    logging.info("is-mapping pipeline completed successfully!")
