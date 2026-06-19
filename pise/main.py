import os
import sys
import argparse
import logging
import concurrent.futures
from pise import filter_reads

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

def get_all_is_elements(primers_file):
    is_elements = set()
    from Bio import SeqIO
    try:
        with open(primers_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split(':')
                if len(parts) >= 2: is_elements.add(parts[1])
    except Exception as e:
        logging.error(f"Error reading primers database file: {e}")
    return sorted(list(is_elements))

def find_matching_is_element(is_name, db_elements):
    for db_element in db_elements:
        if db_element.split('_')[0] == is_name: return db_element
    return None

def discover_samples(input_dir):
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

def str2bool(v):
    if isinstance(v, bool): return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'): return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'): return False
    else: raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_project_root = os.path.dirname(package_dir)
    default_config_path = os.path.join(pkg_project_root, "config", "config.yaml")

    parser = argparse.ArgumentParser(description="pise: Pipeline for IS-Seq sequencing data analysis.")
    parser.add_argument("forward_reads", help="Path to forward reads FASTQ file (_1.fastq.gz) or raw reads directory")
    parser.add_argument("reverse_reads", nargs='?', default=None, help="Path to reverse reads FASTQ file")
    parser.add_argument("-c", "--config", default=None, help="Path to pipeline configuration YAML file")
    parser.add_argument("-t", "--threads", type=int, default=None, help="Number of threads/processes")
    parser.add_argument("--index_i5", type=str2bool, default=None, help="Enable 8 bp i5 index prefix sequence offset")
    parser.add_argument("--qc", type=str2bool, default=None, help="Run QC and demultiplexing step (default: True)")

    args = parser.parse_args()
    config_path = os.path.abspath(args.config if args.config else default_config_path)

    if not os.path.exists(config_path):
        print(f"Error: Config file not found at {config_path}", file=sys.stderr)
        sys.exit(1)
    
    try: config = load_config(config_path)
    except Exception as e:
        print(f"Error parsing configuration file: {e}", file=sys.stderr)
        sys.exit(1)

    target_is_element = config.get("target_is_element", "IS1R_IS1")
    
    primers_file = config.get("primers_file", "config/primers.fasta")
    if primers_file and not os.path.isabs(primers_file):
        config_dir = os.path.dirname(config_path)
        candidate_path = os.path.abspath(os.path.join(os.path.dirname(config_dir), primers_file))
        if os.path.exists(candidate_path): primers_file = candidate_path
        else:
            candidate_path = os.path.abspath(os.path.join(pkg_project_root, primers_file))
            primers_file = candidate_path if os.path.exists(candidate_path) else os.path.abspath(primers_file)

    output_dir = os.path.abspath(config.get("output_dir", "results"))
    os.makedirs(output_dir, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(os.path.join(output_dir, "pise.log"), mode='w')]
    )

    preprocess_cfg = config.get("preprocessing", {})
    min_cov = preprocess_cfg.get("min_cov", 0.98)
    min_identity = preprocess_cfg.get("min_identity", 0.9)
    min_len = preprocess_cfg.get("min_len", 20)
    index_i5 = args.index_i5 if args.index_i5 is not None else preprocess_cfg.get("index_i5", True)
    
    qc_val = preprocess_cfg.get("qc", True)
    if isinstance(qc_val, str):
        if qc_val.lower() in ('none', 'false', 'no', '0'):
            qc_val = False
        else:
            qc_val = True
    elif not isinstance(qc_val, bool):
        qc_val = bool(qc_val)
    qc = args.qc if args.qc is not None else qc_val

    len_actual_primer = preprocess_cfg.get("len_actual_primer", None)
    threads = args.threads if args.threads is not None else preprocess_cfg.get("threads", 4)

    is_batch_mode = os.path.isdir(args.forward_reads)
    batch_tasks = []

    if is_batch_mode:
        raw_samples = discover_samples(args.forward_reads)
        if not raw_samples:
            logging.error(f"No paired FASTQ files found in {args.forward_reads}")
            sys.exit(1)
        
        db_elements = get_all_is_elements(primers_file)
        for sample_prefix, f1_path, f2_path in raw_samples:
            is_name = sample_prefix.split('-')[-1]
            matched_element = find_matching_is_element(is_name, db_elements)
            if matched_element:
                batch_tasks.append({
                    'sample_id': sample_prefix, 'forward': f1_path, 'reverse': f2_path, 'target_is_element': matched_element
                })
            else:
                logging.warning(f"Skipping {sample_prefix}: IS-name {is_name} not found in primers database.")
        
        if not batch_tasks: sys.exit(1)
    else:
        sample_id = os.path.basename(args.forward_reads).replace("_1.fastq.gz", "").replace(".fastq.gz", "")
        db_elements = get_all_is_elements(primers_file)
        is_name = sample_id.split('-')[-1]
        matched_element = find_matching_is_element(is_name, db_elements)
        chosen_target = matched_element if matched_element else target_is_element
        batch_tasks = [{'sample_id': sample_id, 'forward': args.forward_reads, 'reverse': args.reverse_reads, 'target_is_element': chosen_target}]

    logging.info("Initializing IS-Seq Analysis Pipeline (pise)...")
    preprocess_output_dir = os.path.join(output_dir, "filtered_reads")
    os.makedirs(preprocess_output_dir, exist_ok=True)

    qc_results = {}

    # PHASE 1: Parallel QC
    if qc:
        logging.info(f"--- PHASE 1: Parallel QC (Max Threads: {threads}) ---")
        from pise.qc import run_qc
        
        futures = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for task in batch_tasks:
                cur_sample = task['sample_id']
                cur_forward = task['forward']
                cur_target = task['target_is_element']
                
                primers = filter_reads.load_primers(primers_file, cur_target)
                if not primers:
                    logging.error(f"No primers found for {cur_target}. Skipping QC for {cur_sample}.")
                    continue
                
                logging.info(f"Submitting QC task for {cur_sample}")
                f = executor.submit(run_qc, cur_forward, cur_sample, preprocess_output_dir, primers_file, primers, index_i5)
                futures[f] = task
                
        for future in concurrent.futures.as_completed(futures):
            task = futures[future]
            cur_sample = task['sample_id']
            try:
                fw_total, fw_polyn, demux_files = future.result()
                logging.info(f"  QC Finished for {cur_sample}: Total = {fw_total}, Poly-N = {fw_polyn}. Generated {len(demux_files)} files.")
                qc_results[cur_sample] = {
                    'QC_Total_Forward_Reads': fw_total,
                    'QC_PolyN_Forward_Reads': fw_polyn,
                    'files_to_filter': demux_files
                }
            except Exception as exc:
                logging.error(f"QC generated an exception for {cur_sample}: {exc}")

    # PHASE 2: Filtering
    logging.info("--- PHASE 2: Filtering (Internally Parallelized) ---")
    all_runs_stats = []

    for task in batch_tasks:
        cur_sample = task['sample_id']
        cur_forward = task['forward']
        cur_target = task['target_is_element']

        primers = filter_reads.load_primers(primers_file, cur_target)
        if not primers: continue

        qc_stats = {}
        files_to_filter = [cur_forward]
        
        if qc and cur_sample in qc_results:
            qc_data = qc_results[cur_sample]
            qc_stats['QC_Total_Forward_Reads'] = qc_data['QC_Total_Forward_Reads']
            qc_stats['QC_PolyN_Forward_Reads'] = qc_data['QC_PolyN_Forward_Reads']
            files_to_filter = qc_data['files_to_filter']
            
        logging.info(f"--- Filtering sample '{cur_sample}' (Target IS-element: {cur_target}) ---")

        for f_path in files_to_filter:
            if "index-unknown" in os.path.basename(f_path):
                logging.info(f"Skipping filtering for unknown index file: {os.path.basename(f_path)}")
                continue

            # Determine output paths
            base_name = os.path.basename(f_path).replace(".fastq.gz", "")
            if base_name.endswith("_non-polyN_1"):
                base_name = base_name.replace("_non-polyN_1", "_1")
                
            out_full = os.path.join(preprocess_output_dir, f"{base_name}_full.fastq.gz")
            out_bin = os.path.join(preprocess_output_dir, f"{base_name}_bin.fastq.gz")
            out_partial = os.path.join(preprocess_output_dir, f"{base_name}_partial.fastq.gz") if len_actual_primer else None
            
            fw_result = filter_reads.process_forward_reads(
                f_path, primers, min_cov, min_identity, min_len, index_i5, len_actual_primer,
                out_full, out_partial, out_bin, threads=threads
            )
            
            binned_reads = fw_result['total_reads'] - (fw_result['passed_c2_full'] + fw_result['passed_c2_partial'])
            
            run_stats = {
                'File': os.path.basename(f_path),
                'Target_IS_Element': cur_target,
            }
            if qc_stats: run_stats.update(qc_stats)
            
            run_stats.update({
                'Total_Forward_Reads': fw_result['total_reads'],
                'Full_Passed_C1': fw_result['passed_c1_full'],
                'Full_Passed_C2': fw_result['passed_c2_full'],
                'Partial_Passed_C1': fw_result['passed_c1_partial'],
                'Partial_Passed_C2': fw_result['passed_c2_partial'],
                'Binned_Forward_Reads': binned_reads,
            })
            all_runs_stats.append(run_stats)

    pise_summary_path = os.path.join(output_dir, "pise_summary.tsv")
    try:
        os.makedirs(os.path.dirname(pise_summary_path), exist_ok=True)
        with open(pise_summary_path, 'w') as f:
            if all_runs_stats:
                headers = list(all_runs_stats[0].keys())
                f.write("\t".join(headers) + "\n")
                for stats in all_runs_stats:
                    f.write("\t".join([str(stats.get(h, '')) for h in headers]) + "\n")
        logging.info(f"Summary statistics saved to TSV: {pise_summary_path}")
    except Exception as e:
        logging.error(f"Failed to write summary statistics TSV: {e}")

    logging.info("[INFO] To extract reverse reads matching these filtered forward reads, use pise/extract_pairs.py")
    logging.info("Pipeline execution completed successfully!")

if __name__ == "__main__":
    main()
