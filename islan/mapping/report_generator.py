import os
import logging
import csv

try:
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pysam
    from Bio import SeqIO
except ImportError as e:
    logging.warning(f"Required visualization libraries missing: {e}. Falling back to standard reporting.")
    pd = None
    np = None
    go = None
    pysam = None

from islan.constants import (
    MAX_ALIGNMENT_READS,
    GC_WINDOW_SIZE,
    FLANK_PADDING,
    PLOT_HEIGHT,
    READ_PADDING,
    ZOOM_THRESHOLD,
    EXTENSION_PADDING,
    BASE_COLORS,
    HEAD_READ_REVERSE_COLOR,
    HEAD_READ_FORWARD_COLOR,
    TAIL_READ_REVERSE_COLOR,
    TAIL_READ_FORWARD_COLOR,
    HEAD_COV_COLOR,
    TAIL_COV_COLOR,
    GC_LINE_COLOR,
    GENE_FORWARD_COLOR,
    GENE_REVERSE_COLOR,
    MAX_TSD_OVERLAP
)

def parse_bed_cov(filepath):
    cov = []
    if not filepath or not os.path.exists(filepath):
        return cov
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                cov.append((parts[0], int(parts[1]), int(parts[2]), float(parts[3])))
    return cov

def get_base_depths(cov, chrom, start, end):
    depth_map = {}
    for c_chrom, c_start, c_end, c_depth in cov:
        if c_chrom != chrom:
            continue
        if c_start >= end:
            break
        if c_end <= start:
            continue
        overlap_start = max(start, c_start)
        overlap_end = min(end, c_end)
        for pos in range(overlap_start, overlap_end):
            depth_map[pos] = c_depth
    return depth_map

def load_reference_seq(ref_file):
    ref_seqs = {}
    if not ref_file or not os.path.exists(ref_file):
        return ref_seqs
    try:
        if ref_file.endswith('.gbk') or ref_file.endswith('.gb'):
            for record in SeqIO.parse(ref_file, "genbank"):
                ref_seqs[record.id] = str(record.seq)
                ref_seqs[record.name] = str(record.seq)
        else:
            for record in SeqIO.parse(ref_file, "fasta"):
                ref_seqs[record.id] = str(record.seq)
    except Exception as e:
        logging.warning(f"Error loading reference sequence from {ref_file}: {e}")
    return ref_seqs

def get_seq_by_chrom(ref_seqs, chrom):
    if not ref_seqs:
        return None
    if chrom in ref_seqs:
        return ref_seqs[chrom]
    chrom_no_ver = chrom.split('.')[0]
    if chrom_no_ver in ref_seqs:
        return ref_seqs[chrom_no_ver]
    for key in ref_seqs:
        if key.startswith(chrom) or chrom.startswith(key) or key.split('.')[0] == chrom_no_ver:
            return ref_seqs[key]
    return None

def load_genbank_features(ref_file):
    features_dict = {}
    if not ref_file or not os.path.exists(ref_file):
        return features_dict
    if not (ref_file.endswith('.gbk') or ref_file.endswith('.gb')):
        return features_dict
    try:
        for record in SeqIO.parse(ref_file, "genbank"):
            features_list = []
            for feat in record.features:
                if feat.type in ['CDS', 'gene', 'tRNA', 'rRNA']:
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    strand = feat.location.strand
                    locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                    product = feat.qualifiers.get("product", [""])[0]
                    gene_name = feat.qualifiers.get("gene", [""])[0]
                    
                    features_list.append({
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "locus_tag": locus_tag,
                        "product": product,
                        "gene": gene_name,
                        "type": feat.type
                    })
            features_dict[record.id] = features_list
            features_dict[record.id.split('.')[0]] = features_list
    except Exception as e:
        logging.warning(f"Error loading GenBank features from {ref_file}: {e}")
    return features_dict

def get_gc_content_profile(ref_seq, start, end, window_size=GC_WINDOW_SIZE):
    gc_x = list(range(start, end, 10))
    if not gc_x:
        return [], []
    if gc_x[-1] != end - 1:
        gc_x.append(end - 1)
        
    gc_y = []
    half = window_size // 2
    for pos in gc_x:
        sub_seq = ref_seq[max(0, pos - half):min(len(ref_seq), pos + half + 1)]
        if not sub_seq:
            gc_y.append(0.0)
            continue
        gc_count = sum(1 for char in sub_seq if char in 'GCgc')
        gc_pct = (gc_count / len(sub_seq)) * 100
        gc_y.append(gc_pct)
        
    return gc_x, gc_y

def parse_alignments_for_flanks(bam_path, chrom, start, end, ref_seq, flank_label):
    if not pysam or not bam_path or not os.path.exists(bam_path):
        return []
        
    reads_list = []
    try:
        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            bam_chrom = chrom
            if chrom not in samfile.references:
                chrom_no_ver = chrom.split('.')[0]
                for r in samfile.references:
                    if r.split('.')[0] == chrom_no_ver:
                        bam_chrom = r
                        break
            
            if bam_chrom not in samfile.references:
                return []
                
            for read in samfile.fetch(bam_chrom, start, end):
                if read.is_unmapped:
                    continue
                qname = read.query_name
                r_start = read.reference_start
                r_end = read.reference_end
                is_rev = read.is_reverse
                seq = read.query_sequence
                aligned_pairs = read.get_aligned_pairs(with_seq=True)
                
                snps = []
                if seq and ref_seq:
                    for q_pos, r_pos, ref_char in aligned_pairs:
                        if q_pos is not None and r_pos is not None:
                            if start <= r_pos < end:
                                r_base = ref_seq[r_pos]
                                q_base = seq[q_pos]
                                if r_base.upper() != q_base.upper() and r_base.upper() in 'ACGT' and q_base.upper() in 'ACGT':
                                    snps.append({
                                        "pos": r_pos,
                                        "ref": r_base.upper(),
                                        "read": q_base.upper()
                                    })
                                    
                reads_list.append({
                    "qname": qname,
                    "start": r_start,
                    "end": r_end,
                    "is_reverse": is_rev,
                    "snps": snps,
                    "flank": flank_label
                })
    except Exception as e:
        logging.warning(f"Error parsing BAM file {bam_path} in range {start}-{end}: {e}")
        
    return reads_list

def generate_combined_alignment_plotly(left_bam, right_bam, chrom, l_start, l_end, r_start, r_end, ref_seq, features, left_cov, right_cov, cutoff, title, l_label="head", r_label="tail"):
    # GC Content extends outside flanking region
    min_flank = min(l_start, r_start)
    max_flank = max(l_end, r_end)
    start_ext = max(0, min_flank - EXTENSION_PADDING)
    end_ext = max_flank + EXTENSION_PADDING
    if ref_seq:
        end_ext = min(len(ref_seq), end_ext)
        
    # 1. GC content profile across extended range
    gc_x, gc_y = [], []
    if ref_seq:
        gc_x, gc_y = get_gc_content_profile(ref_seq, start_ext, end_ext, window_size=GC_WINDOW_SIZE)
        
    # 2. Left flank coverage profile
    l_depth_map = get_base_depths(left_cov, chrom, l_start, l_end)
    l_coords = sorted(l_depth_map.keys())
    if l_coords:
        l_x = [l_coords[0] - 1] + l_coords + [l_coords[-1] + 1]
        l_y = [0.0] + [l_depth_map[x] for x in l_coords] + [0.0]
    else:
        l_x, l_y = [], []
        
    # 3. Right flank coverage profile
    r_depth_map = get_base_depths(right_cov, chrom, r_start, r_end)
    r_coords = sorted(r_depth_map.keys())
    if r_coords:
        r_x = [r_coords[0] - 1] + r_coords + [r_coords[-1] + 1]
        r_y = [0.0] + [r_depth_map[x] for x in r_coords] + [0.0]
    else:
        r_x, r_y = [], []

    # 4. Mapped reads from left and right BAMs labeled by biological identity (head vs tail)
    left_reads = parse_alignments_for_flanks(left_bam, chrom, l_start, l_end, ref_seq, l_label)
    right_reads = parse_alignments_for_flanks(right_bam, chrom, r_start, r_end, ref_seq, r_label)
    all_reads = left_reads + right_reads
    
    # Group by query name
    templates = {}
    for r in all_reads:
        qname = r["qname"]
        if qname not in templates:
            templates[qname] = []
        templates[qname].append(r)
        
    template_list = []
    for qname, reads in templates.items():
        min_start = min(r["start"] for r in reads)
        max_end = max(r["end"] for r in reads)
        template_list.append({
            "query_name": qname,
            "reads": reads,
            "min_start": min_start,
            "max_end": max_end
        })
        
    template_list.sort(key=lambda t: t["min_start"])
    
    # Downsample templates to MAX_ALIGNMENT_READS (100)
    max_templates = MAX_ALIGNMENT_READS
    if len(template_list) > max_templates:
        step = len(template_list) / max_templates
        downsampled = []
        for i in range(max_templates):
            idx = int(i * step)
            if idx < len(template_list):
                downsampled.append(template_list[idx])
        template_list = downsampled
        
    # Stack templates
    stacked_templates = []
    row_ends = []
    padding = READ_PADDING
    
    for temp in template_list:
        placed = False
        for row_idx, end_coord in enumerate(row_ends):
            if end_coord + padding < temp["min_start"]:
                row_ends[row_idx] = temp["max_end"]
                stacked_templates.append({
                    **temp,
                    "row": row_idx + 1
                })
                placed = True
                break
        if not placed:
            row_ends.append(temp["max_end"])
            stacked_templates.append({
                **temp,
                "row": len(row_ends)
            })
            
    num_rows = max(1, len(row_ends))
    
    # 5. Genes track
    overlapping_features = []
    if features:
        for f in features:
            if f["start"] < end_ext and f["end"] > start_ext:
                overlapping_features.append(f)
                
    # Stack overlapping genes into alternating rows
    gene_rows = []
    for feat in overlapping_features:
        placed = False
        for row_idx, last_end in enumerate(gene_rows):
            if last_end + 100 < feat["start"]:
                gene_rows[row_idx] = feat["end"]
                feat["row"] = row_idx + 1
                placed = True
                break
        if not placed:
            gene_rows.append(feat["end"])
            feat["row"] = len(gene_rows)
            
    num_gene_rows = max(1, len(gene_rows))
    
    # Subplots setup
    has_genes = len(overlapping_features) > 0
    row_heights = [0.35, 0.15, 0.5] if has_genes else [0.45, 0.01, 0.54]
    
    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        row_heights=row_heights,
        specs=[[{"secondary_y": True}], [{"secondary_y": False}], [{"secondary_y": False}]]
    )
    
    # Plot Left Flank Coverage (HEAD in F insertions, TAIL in R insertions)
    if l_x:
        l_cov_name = "HEAD Coverage" if l_label == "head" else "TAIL Coverage"
        l_cov_color = HEAD_COV_COLOR if l_label == "head" else TAIL_COV_COLOR
        fig.add_trace(
            go.Scatter(
                x=l_x, y=l_y,
                name=l_cov_name,
                fill='tozeroy',
                line=dict(color=l_cov_color, width=1.5),
                hovertemplate=f"<b>{l_cov_name}</b><br>Position: %{{x}}<br>Depth: %{{y}}<extra></extra>"
            ),
            row=1, col=1, secondary_y=False
        )
        
    # Plot Right Flank Coverage (TAIL in F insertions, HEAD in R insertions)
    if r_x:
        r_cov_name = "TAIL Coverage" if r_label == "tail" else "HEAD Coverage"
        r_cov_color = TAIL_COV_COLOR if r_label == "tail" else HEAD_COV_COLOR
        fig.add_trace(
            go.Scatter(
                x=r_x, y=r_y,
                name=r_cov_name,
                fill='tozeroy',
                line=dict(color=r_cov_color, width=1.5),
                hovertemplate=f"<b>{r_cov_name}</b><br>Position: %{{x}}<br>Depth: %{{y}}<extra></extra>"
            ),
            row=1, col=1, secondary_y=False
        )
        
    # Plot GC content line (Green, Right Y-Axis)
    if gc_x:
        fig.add_trace(
            go.Scatter(
                x=gc_x, y=gc_y,
                name="GC Content",
                line=dict(color=GC_LINE_COLOR, width=1.2, dash='dot'),
                hovertemplate="<b>GC Content</b><br>Position: %{x}<br>GC %: %{y:.1f}%<extra></extra>"
            ),
            row=1, col=1, secondary_y=True
        )
        
    # Plot Genes
    if has_genes:
        for feat in overlapping_features:
            f_start = max(start_ext, feat["start"])
            f_end = min(end_ext, feat["end"])
            strand_str = "+" if feat["strand"] > 0 else "-"
            color = GENE_FORWARD_COLOR if feat["strand"] > 0 else GENE_REVERSE_COLOR
            label = feat["gene"] or feat["locus_tag"] or "CDS"
            hover_text = f"Gene: {feat['gene']}<br>Locus: {feat['locus_tag']}<br>Product: {feat['product']}<br>Strand: {strand_str}"
            row_y = feat["row"]
            
            fig.add_trace(
                go.Scatter(
                    x=[f_start, f_end],
                    y=[row_y, row_y],
                    mode="lines+text",
                    line=dict(color=color, width=16),
                    text=[label],
                    textposition="middle center",
                    textfont=dict(color="white", size=10, weight="bold"),
                    name=label,
                    hoverinfo="text",
                    hovertext=hover_text,
                    showlegend=False
                ),
                row=2, col=1
            )
            
    # Plot Stacked Reads
    for temp in stacked_templates:
        y = temp["row"]
        reads = temp["reads"]
        
        # Link line if paired
        if len(reads) >= 2:
            r1, r2 = reads[0], reads[1]
            gap_start = min(r1["end"], r2["end"])
            gap_end = max(r1["start"], r2["start"])
            if gap_start < gap_end:
                fig.add_trace(
                    go.Scatter(
                        x=[gap_start, gap_end],
                        y=[y, y],
                        mode="lines",
                        line=dict(color="rgba(160, 160, 160, 0.55)", width=1),
                        showlegend=False,
                        hoverinfo="skip"
                    ),
                    row=3, col=1
                )
                
        # Draw reads
        for read in reads:
            # Color-code based on biological identity: head reads = Blue palette, tail reads = Red palette
            if read["flank"] == "head":
                read_color = HEAD_READ_REVERSE_COLOR if read["is_reverse"] else HEAD_READ_FORWARD_COLOR
            else:
                read_color = TAIL_READ_REVERSE_COLOR if read["is_reverse"] else TAIL_READ_FORWARD_COLOR
                
            fig.add_trace(
                go.Scatter(
                    x=[read["start"], read["end"]],
                    y=[y, y],
                    mode="lines",
                    line=dict(color=read_color, width=8),
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=f"Read ({read['flank'].capitalize()}): {temp['query_name']}<br>Range: {read['start']}-{read['end']}<br>Strand: {'Reverse' if read['is_reverse'] else 'Forward'}"
                ),
                row=3, col=1
            )
            
            # Draw SNPs on this read
            for snp in read["snps"]:
                snp_color = BASE_COLORS.get(snp["read"], "#6b7280")
                fig.add_trace(
                    go.Scatter(
                        x=[snp["pos"]],
                        y=[y],
                        mode="markers+text",
                        marker=dict(color=snp_color, size=6, symbol="square"),
                        text=[snp["read"]],
                        textposition="middle center",
                        textfont=dict(color="white", size=5, family="Courier New, monospace"),
                        showlegend=False,
                        hoverinfo="text",
                        hovertext=f"SNP: Ref {snp['ref']} -> Read {snp['read']}<br>Position: {snp['pos']}"
                    ),
                    row=3, col=1
                )
                
    # 6. Reference Sequence Track
    if ref_seq:
        ref_chars = []
        ref_positions = []
        for p in range(start_ext, end_ext):
            ref_chars.append(ref_seq[p])
            ref_positions.append(p)
            
        fig.add_trace(
            go.Scatter(
                x=ref_positions,
                y=[0] * len(ref_positions),
                mode="text",
                text=ref_chars,
                textfont=dict(size=9, family="Courier New, monospace", color="#424242"),
                name="Ref Seq",
                showlegend=False,
                hoverinfo="skip",
                visible=True if (end_ext - start_ext <= ZOOM_THRESHOLD) else False
            ),
            row=3, col=1
        )
        
    # Configure axes and layout
    fig.update_xaxes(showgrid=True, gridcolor="#f3f4f6", row=1, col=1)
    fig.update_xaxes(showgrid=True, gridcolor="#f3f4f6", row=2, col=1)
    fig.update_xaxes(title_text="Genome Position", showgrid=True, gridcolor="#f3f4f6", range=[start_ext, end_ext], row=3, col=1)
    
    fig.update_yaxes(title_text="Depth", row=1, col=1, secondary_y=False)
    if gc_x:
        fig.update_yaxes(title_text="GC %", range=[0, 100], row=1, col=1, secondary_y=True)
        
    if has_genes:
        fig.update_yaxes(showgrid=False, showticklabels=False, range=[0.2, num_gene_rows + 0.8], row=2, col=1)
    else:
        fig.update_yaxes(showgrid=False, showticklabels=False, row=2, col=1)
        
    fig.update_yaxes(title_text="Stacked Reads", showgrid=False, showticklabels=False, range=[-0.5, num_rows + 0.5], row=3, col=1)
    
    total_height = PLOT_HEIGHT
    
    fig.update_layout(
        title_text=title,
        height=total_height,
        template="plotly_white",
        margin=dict(l=60, r=60, t=50, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    
    return fig

def generate_report(table_file, report_file, reference_file=None, cutoff=6):
    """
    Generates an interactive HTML report using Plotly and Pandas.
    Now includes a combined Coverage/GC subplot, Gene annotations track,
    and stacked alignments (Artemis/IGV-like view) with pair links and SNP highlights.
    Left flank and right flank are shown together on a single plot,
    and the GC content profile extends 5kb on both sides.
    HEAD and TAIL flanking regions are colored distinctly (Blue vs. Red).
    """
    logging.info(f"Generating HTML visualization report: {report_file} (min_y={cutoff})")
    
    if pd is None:
        logging.warning("Required packages for HTML report generation are missing. Skipping.")
        return

    if not os.path.exists(table_file):
        logging.warning(f"Table file {table_file} not found. Skipping report.")
        return

    df = pd.read_csv(table_file, sep='\t')
    if df.empty:
        logging.warning("Table is empty. Skipping report.")
        return

    out_dir = os.path.dirname(os.path.abspath(table_file))
    filename = os.path.basename(table_file)
    
    left_cov_file = None
    right_cov_file = None
    left_bam_file = None
    right_bam_file = None
    
    if '__' in filename:
        sample_prefix, ref_base_with_ext = filename.split('__', 1)
        ref_base = ref_base_with_ext.replace('_table.tsv', '')
        left_cov_file = os.path.join(out_dir, 'tmp', f"{sample_prefix}_left_{ref_base}_cov.bed")
        right_cov_file = os.path.join(out_dir, 'tmp', f"{sample_prefix}_right_{ref_base}_cov.bed")
        left_bam_file = os.path.join(out_dir, f"{sample_prefix}_left_{ref_base}.sorted.bam")
        right_bam_file = os.path.join(out_dir, f"{sample_prefix}_right_{ref_base}.sorted.bam")

    left_cov = parse_bed_cov(left_cov_file)
    right_cov = parse_bed_cov(right_cov_file)
    
    ref_seqs = load_reference_seq(reference_file)
    features_dict = load_genbank_features(reference_file)

    summary = df.groupby('call').size().reset_index(name='count')
    summary_html = summary.to_html(classes='table table-striped', index=False)

    insertions_html = ""
    plot_ids = []
    
    for idx, row in df.iterrows():
        row_num = idx + 1
        chrom = row['contig']
        orient = row.get('orientation', 'F')
        region_name = row.get('region', f"region_{row_num}")
        call_type = row['call']
        
        try:
            l_start, l_end = map(int, str(row['left_pos']).split('-'))
            r_start, r_end = map(int, str(row['right_pos']).split('-'))
        except Exception as e:
            logging.warning(f"Could not parse flank coordinates for row {row_num}: {e}")
            continue

        ref_seq = get_seq_by_chrom(ref_seqs, chrom)
        features = features_dict.get(chrom, [])
        
        # Flank titles and sources depend on predicted orientation
        l_cov_source = left_cov if orient == 'F' else right_cov
        r_cov_source = right_cov if orient == 'F' else left_cov
        
        l_bam_source = left_bam_file if orient == 'F' else right_bam_file
        r_bam_source = right_bam_file if orient == 'F' else left_bam_file
        
        l_label = "head" if orient == "F" else "tail"
        r_label = "tail" if orient == "F" else "head"
        
        title = f"Region {row_num}: {call_type} insertion ({orient}) on {chrom} ({row['left_pos']} to {row['right_pos']})"
        fig = generate_combined_alignment_plotly(
            l_bam_source, r_bam_source, chrom, 
            l_start, l_end, r_start, r_end, 
            ref_seq, features, 
            l_cov_source, r_cov_source, 
            cutoff, title,
            l_label=l_label, r_label=r_label
        )
        
        plot_id = f"plot_r{row_num}_combined"
        plot_ids.append(plot_id)
        plot_html = fig.to_html(full_html=False, include_plotlyjs=False, div_id=plot_id)
        
        insertions_html += f"""
        <div class="card">
            <div class="card-header">
                Region {row_num}: {call_type} insertion on {chrom} ({row['left_pos']} to {row['right_pos']})
            </div>
            <div class="card-body">
                {plot_html}
            </div>
        </div>
        """

    zoom_listener_js = f"""
    <script>
        const plotIds = {plot_ids};
        plotIds.forEach(plotId => {{
            const gd = document.getElementById(plotId);
            if (gd) {{
                gd.on('plotly_relayout', function(eventdata) {{
                    let x0, x1;
                    if (gd.layout && gd.layout.xaxis3 && gd.layout.xaxis3.range) {{
                        x0 = gd.layout.xaxis3.range[0];
                        x1 = gd.layout.xaxis3.range[1];
                    }}
                    if (eventdata['xaxis3.range[0]'] !== undefined) {{
                        x0 = eventdata['xaxis3.range[0]'];
                        x1 = eventdata['xaxis3.range[1]'];
                    }} else if (eventdata['xaxis.range[0]'] !== undefined) {{
                        x0 = eventdata['xaxis.range[0]'];
                        x1 = eventdata['xaxis.range[1]'];
                    }}
                    
                    if (x0 !== undefined && x1 !== undefined) {{
                        const span = x1 - x0;
                        const refSeqIdxs = [];
                        for (let i = 0; i < gd.data.length; i++) {{
                            if (gd.data[i].name === 'Ref Seq') {{
                                refSeqIdxs.push(i);
                            }}
                        }}
                        if (refSeqIdxs.length > 0) {{
                            const show = span <= {ZOOM_THRESHOLD};
                            Plotly.restyle(gd, {{ 'visible': show ? true : false }}, refSeqIdxs);
                        }}
                    }}
                }});
            }}
        }});
    </script>
    """

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IS-Seq Mapping Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 40px;
                background-color: #f3f4f6;
                color: #1f2937;
            }}
            .container {{
                max-width: 1400px;
                margin: auto;
                background: white;
                padding: 35px;
                border-radius: 16px;
                box-shadow: 0 10px 25px rgba(0,0,0,0.05);
            }}
            h1 {{
                color: #111827;
                border-bottom: 2px solid #e5e7eb;
                padding-bottom: 15px;
                font-weight: 700;
                margin-top: 0;
            }}
            h2 {{
                color: #374151;
                margin-top: 35px;
                font-weight: 600;
            }}
            .table {{
                border-collapse: collapse;
                width: 100%;
                margin-top: 15px;
                margin-bottom: 35px;
            }}
            .table th, .table td {{
                border: 1px solid #e5e7eb;
                padding: 14px 16px;
                text-align: left;
            }}
            .table th {{
                background-color: #f9fafb;
                color: #374151;
                font-weight: 600;
            }}
            .table tr:nth-child(even) {{
                background-color: #f9fafb;
            }}
            .card {{
                border: 1px solid #e5e7eb;
                border-radius: 12px;
                margin-bottom: 30px;
                background: white;
                overflow: hidden;
                box-shadow: 0 4px 6px rgba(0,0,0,0.02);
            }}
            .card-header {{
                background-color: #f9fafb;
                padding: 16px 20px;
                font-size: 1.1rem;
                font-weight: 600;
                color: #111827;
                border-bottom: 1px solid #e5e7eb;
            }}
            .card-body {{
                padding: 20px;
            }}
            .warning {{
                color: #b91c1c;
                background: #fef2f2;
                padding: 14px;
                border-left: 4px solid #b91c1c;
                border-radius: 6px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>ISLAN: Insertion Sequence Landscape Analyzer Report</h1>
            <h2>Summary of Detected Insertions</h2>
            {summary_html}
            <div style="font-size: 0.9rem; color: #4b5563; margin-top: -20px; margin-bottom: 25px; background-color: #f9fafb; padding: 12px 16px; border-radius: 8px; border-left: 4px solid #b91c1c;">
                <strong>Note on False Positives (*):</strong> Insertion calls appended with a <strong>*</strong> indicate possible false positives (empty/wild-type loci) where the left and right flanking coverage peaks overlap by more than <strong>{MAX_TSD_OVERLAP} bp</strong>.
            </div>
            <h2>Interactive Flanking Alignment Viewer (IGV/Artemis Style)</h2>
            <p style="color: #6b7280; margin-bottom: 25px;">
                Use the scroll wheel/zoom tools to inspect reads, link lines, and SNPs. 
                Reference sequence characters will dynamically appear when you zoom in (viewport span &le; {ZOOM_THRESHOLD} bp).
            </p>
            {insertions_html}
        </div>
        {zoom_listener_js}
    </body>
    </html>
    """

    with open(report_file, 'w') as f:
        f.write(html_content)
        
    logging.info(f"Report generated successfully.")
