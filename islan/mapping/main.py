import os
import sys
import logging
import argparse
from .mapping import extract_targeted_flanks_pyfastx, map_to_ref_seq, create_bed_files
from .reporting import create_typing_output
from .report_generator import generate_report
from islan.constants import DEFAULT_OUTPUT_DIR, MAPPING_LOG_FILE

def run_is_mapping(args):
    """Entry point for the is-mapping pipeline."""
    # ------------------------------------------------------------------
    # 0. Load optional config file and merge with CLI args.
    #    Priority: CLI arg (not None) > config file > constants default.
    # ------------------------------------------------------------------
    from islan.preprocess import load_config
    from islan.constants import (
        DEFAULT_THREADS, DEFAULT_MIN_DEPTH,
        DEFAULT_MERGING, DEFAULT_MIN_MAPQ, DEFAULT_FLANK_LEN,
        ISElementRegistry,
    )

    cfg_mapping = {}
    if getattr(args, 'config', None):
        try:
            full_cfg = load_config(args.config)
            cfg_mapping = full_cfg.get('is_mapping', {})
        except Exception as exc:
            import sys
            logging.error(f"Failed to load config file '{args.config}': {exc}")
            sys.exit(1)

    def _get(attr, cfg_key, default):
        """Return CLI value if set, else config value, else constant default."""
        cli_val = getattr(args, attr, None)
        if cli_val is not None:
            return cli_val
        if cfg_key in cfg_mapping:
            return cfg_mapping[cfg_key]
        return cfg_mapping.get('cutoff', default)

    # Apply merged values back onto args so the rest of the function is unchanged
    args.threads   = _get('threads',   'threads',   DEFAULT_THREADS)
    args.min_depth = _get('min_depth', 'min_depth', DEFAULT_MIN_DEPTH)
    args.merging   = _get('merging',   'merging',   DEFAULT_MERGING)
    args.min_mapq  = _get('min_mapq',  'min_mapq',  DEFAULT_MIN_MAPQ)
    args.flank_len = _get('flank_len', 'flank_len', DEFAULT_FLANK_LEN)
    # is_name: CLI > config > None
    if not args.is_name:
        args.is_name = cfg_mapping.get('is_name', None) or None

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

    class FlushingFileHandler(logging.FileHandler):
        def emit(self, record):
            super().emit(record)
            self.flush()

    file_handler = FlushingFileHandler(log_file, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    root_logger.addHandler(file_handler)

    sample_prefix = "sample"  # Determine from input

    left_flanking = None
    right_flanking = None

    if args.forward_only:
        if not args.head or not args.tail or not args.filtered_forward:
            logging.error("Forward-only mode requires --head, --tail, and --filtered_forward")
            sys.exit(1)
    else:
        if not args.head or not args.tail or not args.filtered_forward or not args.filtered_reverse:
            logging.error("Targeted mapping requires --head, --tail, --filtered_forward, and --filtered_reverse")
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

    logging.info("Flanking reads extracted. Proceeding to reference mapping.")

    # BWA requires FASTA format. If user provided a GenBank file, convert it to FASTA.
    ref_fasta = args.reference
    if ref_fasta.endswith('.gbk') or ref_fasta.endswith('.gb') or ref_fasta.endswith('.gbff'):
        import Bio.SeqIO
        ref_base_name = os.path.basename(args.reference).rsplit('.', 1)[0]
        ref_fasta = os.path.join(tmp_dir, ref_base_name + '.fasta')
        if not os.path.exists(ref_fasta) or os.path.getsize(ref_fasta) == 0:
            logging.info(f"Converting GenBank to FASTA: {ref_fasta}")
            Bio.SeqIO.convert(args.reference, "genbank", ref_fasta, "fasta")

    # ------------------------------------------------------------------
    # Locate targets.fasta (registry + known-IS scan source).
    # Priority: config targets_file > package config/targets.fasta > primers.fasta
    # ------------------------------------------------------------------
    _pkg_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    cfg_targets = cfg_mapping.get('targets_file', None) or (
        load_config(args.config).get('targets_file', None) if getattr(args, 'config', None) else None
    )

    if cfg_targets:
        # Resolve relative to config file directory if not absolute
        if not os.path.isabs(cfg_targets) and getattr(args, 'config', None):
            cfg_targets = os.path.join(os.path.dirname(os.path.abspath(args.config)), cfg_targets)
        targets_fasta = cfg_targets if os.path.exists(cfg_targets) else None
    else:
        targets_fasta = None

    if not targets_fasta:
        targets_fasta = os.path.join(_pkg_root, "config", "targets.fasta")
    if not os.path.exists(targets_fasta):
        logging.warning(f"Targets file not found. Known IS positions on reference will not be detected.")
        targets_fasta = None



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
        min_depth=args.min_depth,
        merging=args.merging
    )

    ref_base = os.path.basename(args.reference).rsplit('.', 1)[0]
    out_table = os.path.join(out_dir, f"{sample_prefix}__{ref_base}_table.tsv")

    known_is, unpaired_hits = create_typing_output(
        left_merged=left_merged,
        right_merged=right_merged,
        left_cov=os.path.join(tmp_dir, f"{sample_prefix}_left_{ref_base}_cov.bed"),
        right_cov=os.path.join(tmp_dir, f"{sample_prefix}_right_{ref_base}_cov.bed"),
        ref_fasta=args.reference,
        out_file=out_table,
        left_bam=left_bam,
        right_bam=right_bam,
        flank_len=args.flank_len,
        targets_fasta=targets_fasta,
        threads=args.threads,
        is_name=args.is_name
    )

    # Generate HTML report
    report_file = os.path.join(out_dir, f"{sample_prefix}__{ref_base}_report.html")
    generate_report(
        out_table, report_file,
        reference_file=args.reference,
        min_depth=args.min_depth,
        known_is=known_is,
        unpaired_hits=unpaired_hits,
        left_bam=left_bam,
        right_bam=right_bam,
    )

    if not getattr(args, 'temp', False):
        import shutil
        logging.info("Cleaning up temporary files...")
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        logging.info(f"Keeping temporary files in {tmp_dir}")

    logging.info("is-mapping pipeline completed successfully!")


