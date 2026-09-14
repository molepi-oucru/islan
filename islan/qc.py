import os
import sys
import gzip
import logging
from collections import Counter, defaultdict
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from islan.utils import load_known_indices

def generate_qc_report(sample_name, output_dir, total_reads, polyn_reads, length_dist, quality_sums, quality_counts, category_counts, expected_category, primer_variations, primer_hamming, primer_substitutions, index_i5, all_short_names=None):
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
        fig.add_trace(go.Scatter(x=positions, y=avg_quals, mode='lines', name="Avg Quality", line=dict(color='#4cae4c', width=2)), row=2, col=1)
        
        # Add background color zones for quality
        fig.add_hrect(y0=0, y1=20, fillcolor="red", opacity=0.2, layer="below", line_width=0, row=2, col=1)
        fig.add_hrect(y0=20, y1=28, fillcolor="#dca34d", opacity=0.3, layer="below", line_width=0, row=2, col=1)
        fig.add_hrect(y0=28, y1=41, fillcolor="green", opacity=0.2, layer="below", line_width=0, row=2, col=1)
        
        fig.update_yaxes(title_text="Phred Score", range=[0, 40], row=2, col=1)
        fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
    
    # 3. Index Frequencies
    if index_i5 and category_counts:
        if not all_short_names:
            all_short_names = ["IS1R", "ISAeme19", "ISKox3", "ISKpn26"]
        categories = all_short_names + ["undetermined"]
        counts = [category_counts[c] for c in categories]
        
        colors = []
        for c in categories:
            if c == "undetermined":
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
    
    valid_reads = total_reads - polyn_reads
    polyn_pct = (polyn_reads / total_reads * 100) if total_reads > 0 else 0
    valid_pct = (valid_reads / total_reads * 100) if total_reads > 0 else 0
    
    html_content = f"""
    <html>
    <head>
        <title>QC Report - {sample_name}</title>
        <style>
            body {{ font-family: sans-serif; margin: 40px; }}
            .plot-container {{ width: 70%; max-width: 1000px; margin-bottom: 30px; }}
            .heatmap-container {{ width: 315px; height: 315px; margin-bottom: 30px; }}
        </style>
    </head>
    <body>
    <h1>QC Report: {sample_name} (Forward Reads)</h1>
    <ul>
        <li><b>Total Reads:</b> {total_reads}</li>
        <li><b>Poly-N Reads (Discarded):</b> {polyn_reads} ({polyn_pct:.2f}%)</li>
        <li><b>Valid Reads:</b> {valid_reads} ({valid_pct:.2f}%)</li>
    </ul>
    <div class="plot-container">
        {fig.to_html(full_html=False, include_plotlyjs='cdn')}
    </div>
    """
    
    if index_i5:
        html_content += "<h2>Primer Region Quality Control</h2>"
        html_content += "<p><em>Note: This is a fast heuristic classification (reads are officially classified during filtering). The analysis skips the first 8 bases (the index). Only reads matching the targeted/expected index are analyzed for primer variations.</em></p>"
        
        for p_type in ['P_UP', 'P_DOWN']:
            if not primer_hamming[p_type]:
                continue
                
            html_content += f"<h3>{p_type} Primer Analysis</h3>"
            
            # 1. Hamming Distance Distribution
            hamming_counts = primer_hamming[p_type]
            max_h = max(hamming_counts.keys()) if hamming_counts else 0
            x_ham = list(range(max_h + 1))
            y_ham = [hamming_counts.get(i, 0) for i in x_ham]
            
            fig_ham = go.Figure(data=[go.Bar(x=x_ham, y=y_ham)])
            fig_ham.update_layout(title=f"{p_type} Hamming Distance Distribution", xaxis_title="Mismatches", yaxis_title="Reads", height=400)
            html_content += f'<div class="plot-container">{fig_ham.to_html(full_html=False, include_plotlyjs=False)}</div>'
            
            # 2. Base-by-Base Variation
            variations = primer_variations[p_type]
            positions = sorted(variations.keys())
            match_pcts = []
            mismatch_pcts = []
            
            for pos in positions:
                counts = variations[pos]
                total = sum(counts.values())
                match_c = counts.get("Match", 0)
                match_pcts.append((match_c / total * 100) if total > 0 else 0)
                mismatch_pcts.append(100 - match_pcts[-1])
                
            x_pos = [p + 1 for p in positions]
            
            fig_var = go.Figure(data=[
                go.Bar(name='Match %', x=x_pos, y=match_pcts, marker_color='green'),
                go.Bar(name='Mismatch %', x=x_pos, y=mismatch_pcts, marker_color='red')
            ])
            fig_var.update_layout(barmode='stack', title=f"{p_type} Base-by-Base Match Percentage", xaxis_title="Position (1-indexed)", yaxis_title="Percentage", height=400)
            html_content += f'<div class="plot-container">{fig_var.to_html(full_html=False, include_plotlyjs=False)}</div>'
            
            # 3. Nucleotide Substitution Bias
            subs = primer_substitutions[p_type]
            if subs:
                bases = ['A', 'T', 'G', 'C']
                total_subs = sum(subs.values())
                z_data = []
                for read_base in bases: # Y-axis
                    row = []
                    for ref_base in bases: # X-axis
                        if ref_base == read_base:
                            row.append(None)
                        else:
                            key = f"{ref_base}->{read_base}"
                            count = subs.get(key, 0)
                            pct = (count / total_subs * 100) if total_subs > 0 else 0
                            row.append(pct)
                    z_data.append(row)
                    
                fig_sub = go.Figure(data=go.Heatmap(
                    z=z_data,
                    x=bases,
                    y=bases,
                    hoverongaps=False,
                    colorscale='Viridis_r',
                    zmin=0,
                    zmax=100,
                    colorbar=dict(title='%')
                ))
                fig_sub.update_layout(
                    title=f"{p_type} Nucleotide Substitution Bias (%)",
                    xaxis_title="Reference Base (Primer)",
                    yaxis_title="Read Base (Sequenced)",
                    width=450,
                    height=450,
                    yaxis=dict(scaleanchor="x", scaleratio=1)
                )
                html_content += f'<div class="plot-container">{fig_sub.to_html(full_html=False, include_plotlyjs=False)}</div>'
                
    html_content += "</body></html>"
    
    with open(report_path, "w") as f:
        f.write(html_content)

def get_or_create_handle(category, sample_name, output_dir, file_handles, output_files):
    if category not in file_handles:
        if category == "non-polyN":
            out_name = f"{sample_name}_non-polyN_1.fastq.gz"
        elif category == "undetermined":
            out_name = f"{sample_name}_index-undetermined_1.fastq.gz"
        else:
            out_name = f"{sample_name}_index-{category.split('_')[0]}_1.fastq.gz"
        out_path = os.path.join(output_dir, out_name)
        output_files.append(out_path)
        file_handles[category] = gzip.open(out_path, "wt")
    return file_handles[category]

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def run_qc(fastq_path, sample_name, output_dir, primers_fasta_path, primers=None, index_i5=True, i5_mismatch=2, all_short_names=None):
    """
    Runs QC and demultiplexes forward reads by index. Poly-N reads are dropped.
    Returns total_reads, polyn_reads, category_counts, and created output fastq files.

    all_short_names: list of IS short names from ISElementRegistry (e.g. ['IS1R', 'ISKpn26']).
                     If None, falls back to the four hard-coded names for backwards compatibility.
    """
    # Ensure we always have a list to work with
    if not all_short_names:
        all_short_names = ["IS1R", "ISAeme19", "ISKox3", "ISKpn26"]
    # E.g. sample_name: "278-IS1R" -> expected_category: "IS1R"
    expected_category = sample_name.split('-')[-1]
    
    total_reads = 0
    polyn_reads = 0
    
    length_dist = Counter()
    quality_sums = defaultdict(int)
    quality_counts = defaultdict(int)
    
    index_len = 8 if index_i5 else 0
    index_counts = Counter()
    category_counts = Counter({n: 0 for n in all_short_names + ['undetermined']})
    
    # Load the known indices from primers.fasta
    known_indices = load_known_indices(primers_fasta_path) if index_i5 else {}
    
    primer_variations = {
        'P_UP': defaultdict(lambda: Counter()),
        'P_DOWN': defaultdict(lambda: Counter())
    }
    primer_hamming = {
        'P_UP': Counter(),
        'P_DOWN': Counter()
    }
    primer_substitutions = {
        'P_UP': Counter(),
        'P_DOWN': Counter()
    }
    
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
                # Downsample quality profiling (1 in 10 reads)
                if total_reads % 10 == 0:
                    for i, q in enumerate(qual_str):
                        quality_sums[i] += ord(q) - 33
                        quality_counts[i] += 1
                
                category = "non-polyN"
                category_prefix = "non-polyN"
                if index_i5 and len(seq) >= index_len:
                    idx_seq = seq[:index_len]
                    index_counts[idx_seq] += 1
                    
                    # Find closest known index sequence within mismatch tolerance
                    best_dist = i5_mismatch + 1
                    best_category = None
                    best_prefix = None
                    
                    for is_element, kn_seq in known_indices.items():
                        dist = hamming_distance(idx_seq, kn_seq)
                        if dist < best_dist:
                            best_dist = dist
                            best_category = is_element
                            best_prefix = is_element.split('_')[0]
                            
                    if best_category is not None:
                        category = best_category
                        category_prefix = best_prefix
                    else:
                        category = "undetermined"
                        category_prefix = "undetermined"
                    
                    category_counts[category_prefix] += 1
                    
                    # Base variation calculation - only for the true expected category prefix
                    if primers and category_prefix == expected_category:
                        best_p_type = None
                        best_p_mismatches = float('inf')
                        best_sub_seq = None
                        
                        # Calculate hamming distance to find the best matching primer
                        for p_type in ['P_UP', 'P_DOWN']:
                            if p_type not in primers: continue
                            p_seq = primers[p_type]
                            expected_len = len(p_seq)
                            sub_seq = seq[index_len : index_len + expected_len]
                            if len(sub_seq) < expected_len:
                                continue
                            
                            mismatches = sum(1 for i in range(expected_len) if sub_seq[i] != p_seq[i])
                            if mismatches < best_p_mismatches:
                                best_p_mismatches = mismatches
                                best_p_type = p_type
                                best_sub_seq = sub_seq
                                
                        if best_p_type is not None:
                            p_seq = primers[best_p_type]
                            primer_hamming[best_p_type][best_p_mismatches] += 1
                            for pos in range(len(p_seq)):
                                ref_char = p_seq[pos]
                                read_char = best_sub_seq[pos]
                                if read_char == ref_char:
                                    primer_variations[best_p_type][pos]["Match"] += 1
                                else:
                                    primer_variations[best_p_type][pos][f"Mismatch ({read_char})"] += 1
                                    primer_substitutions[best_p_type][f"{ref_char}->{read_char}"] += 1
                                    
                # Write to the assigned demux handle
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
            
    # Generate report in qc_reports directory
    qc_dir = os.path.join(os.path.dirname(output_dir), "qc_reports")
    generate_qc_report(sample_name, qc_dir, total_reads, polyn_reads,
                        length_dist, quality_sums, quality_counts, category_counts,
                        expected_category, primer_variations, primer_hamming, primer_substitutions, index_i5,
                        all_short_names=all_short_names)
                       
    return total_reads, polyn_reads, category_counts, output_files
