import os
import sys
import gzip
import logging
from collections import Counter, defaultdict
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import SeqIO

def load_known_indices(primers_path):
    known_indices = {}
    try:
        with open(primers_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id.endswith(":P_UP"):
                    parts = record.id.split(':')
                    if len(parts) >= 2:
                        is_element = parts[1]
                        known_indices[is_element] = str(record.seq[:8])
    except Exception as e:
        logging.warning(f"Could not load known indices from {primers_path}: {e}")
    return known_indices

def generate_qc_report(sample_name, output_dir, total_reads, polyn_reads, length_dist, quality_sums, quality_counts, category_counts, expected_category, primer_variations, index_i5):
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, f"{sample_name}_Forward_QC_Report.html")
    
    fig = make_subplots(
        rows=3 if index_i5 else 2, cols=1,
        subplot_titles=(
            f"Read Length Distribution",
            f"Average Per-Base Quality",
            f"Index Sequences Frequency" if index_i5 else ""
        ),
        vertical_spacing=0.1
    )
    
    # 1. Length Distribution
    lengths = sorted(length_dist.keys())
    counts = [length_dist[l] for l in lengths]
    if lengths:
        fig.add_trace(go.Bar(x=lengths, y=counts, name="Length Dist"), row=1, col=1)
    
    # 2. Quality Scores
    positions = sorted(quality_sums.keys())
    avg_quals = [quality_sums[p] / quality_counts[p] for p in positions if quality_counts[p] > 0]
    if positions:
        fig.add_trace(go.Scatter(x=positions, y=avg_quals, mode='lines+markers', name="Avg Quality"), row=2, col=1)
    
    # 3. Index Frequencies
    if index_i5 and category_counts:
        categories = ["IS1R", "ISAeme19", "ISKox3", "ISKpn26", "unknown"]
        counts = [category_counts[c] for c in categories]
        
        colors = []
        for c in categories:
            if c == "unknown":
                colors.append("red")
            elif c == expected_category:
                colors.append("green")
            else:
                colors.append("yellow")
        
        fig.add_trace(go.Bar(x=categories, y=counts, marker_color=colors, name="Index Freq"), row=3, col=1)

    fig.update_layout(
        title_text=f"QC Report: {sample_name} (Forward Reads)<br>Total: {total_reads} | Poly-N: {polyn_reads}",
        height=900 if index_i5 else 600,
        showlegend=False
    )
    
    html_content = f"""
    <html>
    <head><title>QC Report - {sample_name}</title></head>
    <body style="font-family: sans-serif; margin: 40px;">
    <h1>QC Report: {sample_name} (Forward Reads)</h1>
    <ul>
        <li><b>Total Reads:</b> {total_reads}</li>
        <li><b>Poly-N Reads (Discarded):</b> {polyn_reads}</li>
        <li><b>Valid Reads:</b> {total_reads - polyn_reads}</li>
    </ul>
    """
    
    html_content += fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    if index_i5:
        html_content += "<h2>Base-by-Base Variation (Primer Region)</h2>"
        html_content += "<p>Matches to the expected primer sequence after the index.</p>"
        for p_type, variations in primer_variations.items():
            if not variations: continue
            html_content += f"<h3>{p_type} Sequence</h3>"
            html_content += "<table border='1' cellpadding='5' style='border-collapse: collapse;'>"
            html_content += "<tr><th>Position</th><th>Match %</th><th>Mismatch Details</th></tr>"
            
            for pos in sorted(variations.keys()):
                counts = variations[pos]
                total = sum(counts.values())
                match_count = counts.get("Match", 0)
                match_pct = (match_count / total * 100) if total > 0 else 0
                
                mismatches = [f"{k}: {v}" for k, v in counts.items() if k != "Match"]
                mismatch_str = ", ".join(mismatches) if mismatches else "None"
                
                html_content += f"<tr><td>{pos+1}</td><td>{match_pct:.1f}%</td><td>{mismatch_str}</td></tr>"
            html_content += "</table>"
            
    html_content += "</body></html>"
    
    with open(report_path, "w") as f:
        f.write(html_content)

def get_or_create_handle(category, sample_name, output_dir, file_handles, output_files):
    if category not in file_handles:
        # File name pattern: sample_indexCategory_1.fastq.gz
        if category == "non-polyN":
            out_name = f"{sample_name}_non-polyN_1.fastq.gz"
        else:
            out_name = f"{sample_name}_index-{category.split('_')[0]}_1.fastq.gz"
        out_path = os.path.join(output_dir, out_name)
        output_files.append(out_path)
        file_handles[category] = gzip.open(out_path, "wt")
    return file_handles[category]

def run_qc(fastq_path, sample_name, output_dir, primers_fasta_path, primers=None, index_i5=True):
    """
    Runs QC and demultiplexes forward reads by index. Poly-N reads are dropped.
    Returns a list of created output fastq files for downstream filtering.
    """
    expected_category = sample_name.split('-')[-1]
    
    total_reads = 0
    polyn_reads = 0
    
    length_dist = Counter()
    quality_sums = defaultdict(int)
    quality_counts = defaultdict(int)
    
    index_len = 8 if index_i5 else 0
    index_counts = Counter()
    category_counts = Counter({
        "IS1R": 0,
        "ISAeme19": 0,
        "ISKox3": 0,
        "ISKpn26": 0,
        "unknown": 0
    })
    
    # Load the 4 known indices from primers.fasta
    known_indices = load_known_indices(primers_fasta_path) if index_i5 else {}
    # Cache mapping from index sequence to (category_full_name, category_prefix) to avoid string splits and loops
    known_indices_map = {kn_seq: (is_element, is_element.split('_')[0]) for is_element, kn_seq in known_indices.items()}
    
    # expected_categories determined from sample name
            
    primer_variations = {
        'HEAD': defaultdict(lambda: Counter()),
        'TAIL': defaultdict(lambda: Counter())
    }
    
    # We will output up to 5 demultiplexed files directly in output_dir
    # e.g., output_dir is preprocess_output_dir (results/filtered_reads)
    output_files = []
    file_handles = {}

    try:
        with gzip.open(fastq_path, "rt") as infile:
            while True:
                header = infile.readline()
                if not header: break
                seq = infile.readline().strip().upper()
                spacer = infile.readline()
                qual_str = infile.readline().strip()
                
                total_reads += 1
                if seq and seq[0] == 'N' and seq == 'N' * len(seq):
                    polyn_reads += 1
                    continue
                
                length_dist[len(seq)] += 1
                # Downsample quality profiling (1 in 10 reads) for a 10x speedup in the inner char loop
                if total_reads % 10 == 0:
                    for i, q in enumerate(qual_str):
                        quality_sums[i] += ord(q) - 33
                        quality_counts[i] += 1
                
                category = "non-polyN"
                category_prefix = "non-polyN"
                if index_i5 and len(seq) >= index_len:
                    idx_seq = seq[:index_len]
                    index_counts[idx_seq] += 1
                    
                    # Demultiplex routing using precomputed map
                    if idx_seq in known_indices_map:
                        category, category_prefix = known_indices_map[idx_seq]
                    else:
                        category = "unknown"
                        category_prefix = "unknown"
                    
                    category_counts[category_prefix] += 1
                    
                    # Base variation calculation - only for the true expected category prefix
                    if primers and category_prefix == expected_category:
                        for p_type, p_seq in primers.items():
                            # Match the start of the primer sequence portion of the read
                            if len(seq) >= index_len + 8 and seq[index_len : index_len + 8] == p_seq[:8]:
                                sub_seq = seq[index_len : index_len + len(p_seq)]
                                for pos in range(min(len(sub_seq), len(p_seq))):
                                    match_char = "Match" if sub_seq[pos] == p_seq[pos] else f"Mismatch ({sub_seq[pos]})"
                                    primer_variations[p_type][pos][match_char] += 1
                                    
                # Write natively to the assigned demux handle
                # We skip SeqIO for massive speedup
                handle = get_or_create_handle(category, sample_name, output_dir, file_handles, output_files)
                handle.write(f"{header}{seq}\n+\n{qual_str}\n")
                
    except FileNotFoundError:
        logging.error(f"QC input file not found at: {fastq_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred during QC of '{fastq_path}': {e}")
        sys.exit(1)
    finally:
        for fh in file_handles.values():
            fh.close()
            
    # Generate report in qc_reports directory (one level up from output_dir if output_dir is filtered_reads)
    # Actually, we should output report to a dedicated qc_reports dir
    qc_dir = os.path.join(os.path.dirname(output_dir), "qc_reports")
    generate_qc_report(sample_name, qc_dir, total_reads, polyn_reads, 
                       length_dist, quality_sums, quality_counts, category_counts, expected_category, primer_variations, index_i5)
                       
    # If index_counts has only 1 unique index (or none), and category is "unknown" or known, it's fine.
    # The handles are already closed.
    
    return total_reads, polyn_reads, output_files
