import os
import sys
import logging
import concurrent.futures
from Bio import SeqIO
from islan import filter_reads
from islan.constants import (
    INDEX_LEN,
    ISElementRegistry,
    DEFAULT_MIN_LEN,
    DEFAULT_I5_MISMATCH,
    DEFAULT_MIN_COV,
    DEFAULT_MIN_IDENTITY,
    DEFAULT_QC,
)


def load_config(filepath):
    config = {}
    current_section = None
    with open(filepath, 'r') as f:
        for line in f:
            indent = len(line) - len(line.lstrip())
            stripped = line.strip()
            if not stripped or stripped.startswith('#'): continue
            if ':' not in stripped: continue
            
            key, val = stripped.split(':', 1)
            key = key.strip()
            val = val.strip()
            if '#' in val: val = val.split('#', 1)[0].strip()
            
            if val == '':
                current_section = key
                config[current_section] = {}
                continue
            
            if val.lower() == 'null': parsed_val = None
            elif val.lower() == 'true': parsed_val = True
            elif val.lower() == 'false': parsed_val = False
            elif (val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'")):
                parsed_val = val[1:-1]
            else:
                try:
                    if '.' in val: parsed_val = float(val)
                    else: parsed_val = int(val)
                except ValueError:
                    parsed_val = val
            
            if indent > 0 and current_section: config[current_section][key] = parsed_val
            else:
                current_section = None
                config[key] = parsed_val
    return config


def discover_samples(input_dir):
    """Discover paired FASTQ sample files in a directory."""
    samples = []
    for filename in sorted(os.listdir(input_dir)):
        if filename.endswith("_1.fastq.gz"):
            f1_path = os.path.join(input_dir, filename)
            f2_filename = filename.replace("_1.fastq.gz", "_2.fastq.gz")
            f2_path = os.path.join(input_dir, f2_filename)
            if os.path.exists(f2_path):
                sample_prefix = filename.replace("_1.fastq.gz", "")
                samples.append((sample_prefix, f1_path, f2_path))
    return samples


def run_preprocess(forward_reads, reverse_reads=None, config_path=None, threads=None, index_i5=None, qc=None, i5_mismatch=None, min_len=None, output_dir=None):

    package_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_project_root = os.path.dirname(package_dir)
    default_config_path = os.path.join(pkg_project_root, "config", "config.yaml")

    config_filepath = os.path.abspath(config_path if config_path else default_config_path)

    if not os.path.exists(config_filepath):
        logging.error(f"Config file not found at {config_filepath}")
        sys.exit(1)
    
    try: config = load_config(config_filepath)
    except Exception as e:
        logging.error(f"Error parsing configuration file: {e}")
        sys.exit(1)

    target_is_element = config.get("target_is_element", None)

    targets_file = config.get("targets_file", config.get("primers_file", "config/targets.fasta"))
    if targets_file and not os.path.isabs(targets_file):
        config_dir = os.path.dirname(config_filepath)
        candidate_path = os.path.abspath(os.path.join(os.path.dirname(config_dir), targets_file))
        if os.path.exists(candidate_path):
            targets_file = candidate_path
        else:
            candidate_path = os.path.abspath(os.path.join(pkg_project_root, targets_file))
            targets_file = candidate_path if os.path.exists(candidate_path) else os.path.abspath(targets_file)
    primers_file = targets_file

    # Build the IS element registry from targets.fasta (single source of truth)
    registry = ISElementRegistry(targets_file)
    if not registry.all_short_names:
        logging.error("Could not build IS element registry from targets_file. "
                      "Check that targets.fasta exists and has valid headers.")
        sys.exit(1)
    logging.info(f"IS element registry loaded: {registry.all_full_names}")


    out_dir_cfg = config.get("output_dir", "results")
    final_output_dir = os.path.abspath(output_dir if output_dir is not None else out_dir_cfg)
    os.makedirs(final_output_dir, exist_ok=True)
    
    # Configure logging immediately
    log_file = os.path.join(final_output_dir, "islan.log")
    class FlushingFileHandler(logging.FileHandler):
        def emit(self, record):
            super().emit(record)
            self.flush()

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    file_h = FlushingFileHandler(log_file, mode='w')
    file_h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    root_logger.addHandler(file_h)

    preprocess_cfg = config.get("preprocessing", {})
    min_cov = preprocess_cfg.get("min_cov", DEFAULT_MIN_COV)
    min_identity = preprocess_cfg.get("min_identity", DEFAULT_MIN_IDENTITY)
    min_len_val = min_len if min_len is not None else preprocess_cfg.get("min_len", DEFAULT_MIN_LEN)
    index_i5_val = index_i5 if index_i5 is not None else preprocess_cfg.get("index_i5", True)
    i5_mismatch_val = i5_mismatch if i5_mismatch is not None else preprocess_cfg.get("i5_mismatch", DEFAULT_I5_MISMATCH)

    qc_val = preprocess_cfg.get("qc", DEFAULT_QC)
    if isinstance(qc_val, str):
        if qc_val.lower() in ('none', 'false', 'no', '0'):
            qc_val = False
        else:
            qc_val = True
    elif not isinstance(qc_val, bool):
        qc_val = bool(qc_val)
    run_qc_step = qc if qc is not None else qc_val

    len_actual_primer = preprocess_cfg.get("len_actual_primer", None)
    threads_val = threads if threads is not None else preprocess_cfg.get("threads", 4)

    is_batch_mode = os.path.isdir(forward_reads)
    batch_tasks = []

    if is_batch_mode:
        raw_samples = discover_samples(forward_reads)
        if not raw_samples:
            logging.error(f"No paired FASTQ files found in {forward_reads}")
            sys.exit(1)
        
        for sample_prefix, f1_path, f2_path in raw_samples:
            matched_element, expected_cat = registry.resolve(sample_prefix)
            if matched_element:
                batch_tasks.append({
                    'sample_id': sample_prefix,
                    'forward': f1_path,
                    'reverse': f2_path,
                    'target_is_element': matched_element,
                    'expected_category': expected_cat
                })
            else:
                logging.warning(f"Skipping {sample_prefix}: Expected IS element name not recognized in filename.")
        
        if not batch_tasks: sys.exit(1)
    else:
        sample_id = os.path.basename(forward_reads).replace("_1.fastq.gz", "").replace(".fastq.gz", "")
        matched_element, expected_cat = registry.resolve(sample_id)
        chosen_target = matched_element if matched_element else (target_is_element or registry.all_full_names[0])
        chosen_cat = expected_cat if expected_cat else chosen_target.split('_')[0]
        
        rev_path = reverse_reads
        if not rev_path and forward_reads.endswith("_1.fastq.gz"):
            candidate_rev = forward_reads.replace("_1.fastq.gz", "_2.fastq.gz")
            if os.path.exists(candidate_rev):
                rev_path = candidate_rev
                
        batch_tasks = [{
            'sample_id': sample_id, 
            'forward': forward_reads, 
            'reverse': rev_path, 
            'target_is_element': chosen_target,
            'expected_category': chosen_cat
        }]

    logging.info("Initializing IS-Seq Analysis Pipeline (islan) pre-process...")
    preprocess_output_dir = os.path.join(final_output_dir, "filtered_reads")
    os.makedirs(preprocess_output_dir, exist_ok=True)
    qc_output_dir = os.path.join(final_output_dir, "demux_reads")
    os.makedirs(qc_output_dir, exist_ok=True)

    qc_results = {}

    # PHASE 1: Parallel QC
    if run_qc_step:
        logging.info(f"--- PHASE 1: Parallel QC (Max Threads: {threads_val}) ---")
        from islan.qc import run_qc
        
        futures = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads_val) as executor:
            for task in batch_tasks:
                cur_sample = task['sample_id']
                cur_forward = task['forward']
                cur_target = task['target_is_element']
                
                primers = filter_reads.load_primers(primers_file, cur_target)
                if not primers:
                    logging.error(f"No primers found for {cur_target}. Skipping QC for {cur_sample}.")
                    continue
                
                logging.info(f"Submitting QC task for {cur_sample}")
                f = executor.submit(run_qc, cur_forward, cur_sample, qc_output_dir, primers_file, primers,
                                   index_i5_val, i5_mismatch_val, registry.all_short_names)
                futures[f] = task

        for future in concurrent.futures.as_completed(futures):
            task = futures[future]
            cur_sample = task['sample_id']
            try:
                fw_total, fw_polyn, category_counts, demux_files = future.result()
                polyn_pct = (fw_polyn / fw_total * 100) if fw_total > 0 else 0
                fw_valid = fw_total - fw_polyn
                valid_pct = (fw_valid / fw_total * 100) if fw_total > 0 else 0

                logging.info(f"  QC Finished for {cur_sample}:")
                logging.info(f"    - Total reads: {fw_total}")
                logging.info(f"    - Poly-N reads: {fw_polyn} ({polyn_pct:.2f}%)")
                logging.info(f"    - Valid reads: {fw_valid} ({valid_pct:.2f}%)")

                qc_results[cur_sample] = {
                    'QC_Total_Forward_Reads': fw_total,
                    'QC_PolyN_Forward_Reads': fw_polyn,
                    'category_counts': category_counts,
                    'files_to_filter': demux_files
                }
            except Exception as exc:
                logging.error(f"QC generated an exception for {cur_sample}: {exc}")
                sys.exit(1)

    # PHASE 2: Filtering
    logging.info("--- PHASE 2: Filtering (Internally Parallelized) ---")
    all_runs_stats = []

    for task in batch_tasks:
        cur_sample = task['sample_id']
        cur_forward = task['forward']
        cur_target = task['target_is_element']
        expected_cat = task['expected_category']

        primers = filter_reads.load_primers(primers_file, cur_target)
        if not primers: 
            logging.error(f"No primers found for {cur_target}. Skipping filtering for {cur_sample}.")
            continue

        head_len = len(primers.get('HEAD', ''))
        tail_len = len(primers.get('TAIL', ''))
        max_primer_len = max(head_len, tail_len)
        min_len_forward = INDEX_LEN + max_primer_len + min_len_val

        logging.info(f"Sample: {cur_sample} | Expected IS: {cur_target} | HEAD len: {head_len} | TAIL len: {tail_len} | MIN_LEN_FORWARD: {min_len_forward}")

        status = "Index OK"
        total_non_polyN = 0
        undetermined_count = 0
        qc_total = 0
        qc_polyn = 0
        
        if run_qc_step and cur_sample in qc_results:
            qc_data = qc_results[cur_sample]
            qc_total = qc_data['QC_Total_Forward_Reads']
            qc_polyn = qc_data['QC_PolyN_Forward_Reads']
            cat_counts = qc_data['category_counts']
            total_non_polyN = qc_total - qc_polyn
            undetermined_count = cat_counts.get("undetermined", 0)

            if total_non_polyN > 0:
                expected_pct = (cat_counts.get(expected_cat, 0) / total_non_polyN) * 100
                if expected_pct < 30.0:
                    status = "Index failure"
                else:
                    other_warning = False
                    for other_cat in registry.all_short_names:
                        if other_cat != expected_cat:
                            other_pct = (cat_counts.get(other_cat, 0) / total_non_polyN) * 100
                            if other_pct >= 30.0:
                                other_warning = True
                                break
                    status = "Index warning" if other_warning else "Index OK"
            else:
                status = "Index failure"

        logging.info(f"Index Status for {cur_sample}: {status}")

        if status == "Index failure":
            logging.warning(f"Skipping filtering for {cur_sample} due to Index failure.")
            run_stats = {
                'Sample_Index': cur_sample,
                'IS_element': cur_target,
                'QC_Total_Forward_Reads': qc_total,
                'QC_PolyN_Forward_Reads': qc_polyn,
                'QC_Undetermined_Forward_Reads': undetermined_count,
                'QC_NonPolyN_Forward_Reads': total_non_polyN,
                'Filtering_Status': status
            }
            all_runs_stats.append(run_stats)
            continue

        filter_input_file = cur_forward
        expected_index_seq = None
        if run_qc_step and cur_sample in qc_results:
            expected_demux_file = os.path.join(qc_output_dir, f"{cur_sample}_index-{expected_cat}_1.fastq.gz")
            if os.path.exists(expected_demux_file):
                filter_input_file = expected_demux_file
            else:
                logging.warning(f"Demultiplexed file for {expected_cat} not found at {expected_demux_file}. Using raw forward reads.")
        else:
            from islan.qc import load_known_indices
            known_indices = load_known_indices(primers_file)
            expected_index_seq = known_indices.get(cur_target, None)

        out_filt = os.path.join(preprocess_output_dir, f"{cur_sample}_filtered_1.fastq.gz")
        out_head = os.path.join(preprocess_output_dir, f"{cur_sample}_filtered_1_HEAD.fastq.gz")
        out_tail = os.path.join(preprocess_output_dir, f"{cur_sample}_filtered_1_TAIL.fastq.gz")

        filter_res = filter_reads.process_forward_reads(
            filter_input_file, primers, min_cov, min_identity, min_len_val, INDEX_LEN, min_len_forward,
            out_filt, out_head, out_tail,
            expected_index_seq=expected_index_seq, i5_mismatch=i5_mismatch_val,
            len_actual_primer=len_actual_primer, threads=threads_val
        )

        run_stats = {
            'Sample_Index': cur_sample,
            'IS_element': cur_target,
            'QC_Total_Forward_Reads': qc_total if run_qc_step else filter_res['total_filtered'],
            'QC_PolyN_Forward_Reads': qc_polyn if run_qc_step else 0,
            'QC_Undetermined_Forward_Reads': undetermined_count,
            'QC_NonPolyN_Forward_Reads': total_non_polyN if run_qc_step else filter_res['total_filtered'],
            'Filtering_Status': status
        }
        all_runs_stats.append(run_stats)

    islan_summary_path = os.path.join(final_output_dir, "islan_summary.tsv")
    try:
        with open(islan_summary_path, 'w') as f:
            if all_runs_stats:
                headers = [
                    'Sample_Index', 'IS_element',
                    'QC_Total_Forward_Reads', 'QC_PolyN_Forward_Reads',
                    'QC_Undetermined_Forward_Reads', 'QC_NonPolyN_Forward_Reads',
                    'Filtering_Status'
                ]
                f.write("\t".join(headers) + "\n")
                for stats in all_runs_stats:
                    f.write("\t".join([str(stats.get(h, '')) for h in headers]) + "\n")
        logging.info(f"Summary statistics saved to TSV: {islan_summary_path}")
    except Exception as e:
        logging.error(f"Failed to write summary statistics TSV: {e}")
        sys.exit(1)

    logging.info("Pipeline pre-process execution completed successfully!")
