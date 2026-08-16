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
    MAX_TSD_OVERLAP,
    ALL_REPORT_CALL_CLASSES,
)

# Display-name mapping: maps internal call strings (as written to TSV) to
# human-readable labels shown only in the HTML report.
_DISPLAY_CLASS_NAMES = {
    'Off-Target Amplicon (Noise)': 'Left-Right Imbalance Depth',
}

def _display_call(call_str):
    """Map an internal call string to its HTML display label."""
    return _DISPLAY_CLASS_NAMES.get(call_str, call_str)

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
        if ref_file.endswith('.gbk') or ref_file.endswith('.gb') or ref_file.endswith('.gbff'):
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
    if not (ref_file.endswith('.gbk') or ref_file.endswith('.gb') or ref_file.endswith('.gbff')):
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
        flank_type = reads[0]["flank"] if reads else "head"
        template_list.append({
            "query_name": qname,
            "reads": reads,
            "min_start": min_start,
            "max_end": max_end,
            "flank": flank_type
        })
        
    # Sort templates: Group by flank type (HEAD vs TAIL), then prioritize junction-defining reads (R1)
    def is_junction_template(t):
        for r in t["reads"]:
            if abs(r["start"] - l_start) <= 50 or abs(r["end"] - l_end) <= 50:
                return True
            if abs(r["start"] - r_start) <= 50 or abs(r["end"] - r_end) <= 50:
                return True
        return False

    template_list.sort(key=lambda t: (
        0 if t["flank"] == l_label else 1,
        0 if is_junction_template(t) else 1,
        t["min_start"]
    ))
    
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
    
    # Add vertical highlight for TSD / insertion site overlap
    tsd_candidates = [l_start, l_end, r_start, r_end]
    min_j = min(l_start, r_start)
    max_j = max(l_start, r_start)
    if abs(min_j - max_j) <= 100:
        fig.add_vrect(
            x0=min_j, x1=max_j if max_j != min_j else min_j + 5,
            fillcolor="rgba(251, 191, 36, 0.22)",
            line_color="rgba(217, 119, 6, 0.6)",
            line_width=1,
            line_dash="dot",
            annotation_text="TSD Zone",
            annotation_position="top left",
            row=1, col=1
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
            
    # Plot Stacked Reads (Fragment-Centric Overlay View)
    for temp in stacked_templates:
        y = temp["row"]
        reads = temp["reads"]
        min_start = temp["min_start"]
        max_end = temp["max_end"]
        
        # 1. Draw Full DNA Fragment Span as a thin background line
        fig.add_trace(
            go.Scatter(
                x=[min_start, max_end],
                y=[y, y],
                mode="lines",
                line=dict(color="rgba(156, 163, 175, 0.45)", width=2),
                showlegend=False,
                hoverinfo="skip"
            ),
            row=3, col=1
        )
                
        # 2. Draw reads overlaid on top of the fragment line
        for read in reads:
            # Determine if this read is a junction-defining Read 1 (starts or ends near peak junction)
            is_j_read = (abs(read["start"] - l_start) <= 50 or abs(read["end"] - l_end) <= 50 or 
                          abs(read["start"] - r_start) <= 50 or abs(read["end"] - r_end) <= 50)
            
            # Color-code based on read identity: Junction Read 1 = Solid Dark Color, Flank Read 2 = Light Color
            if read["flank"] == "head":
                read_color = HEAD_READ_FORWARD_COLOR if is_j_read else HEAD_READ_REVERSE_COLOR
            else:
                read_color = TAIL_READ_FORWARD_COLOR if is_j_read else TAIL_READ_REVERSE_COLOR
                
            read_width = 9 if is_j_read else 5
            read_file_str = "_1 (Forward Read)" if is_j_read else "_2 (Reverse Mate)"
            read_type_str = "IS Junction Boundary" if is_j_read else "Extended Flank"
            strand_str = "(+) Forward Strand" if not read["is_reverse"] else "(-) Reverse Strand"
            
            fig.add_trace(
                go.Scatter(
                    x=[read["start"], read["end"]],
                    y=[y, y],
                    mode="lines",
                    line=dict(color=read_color, width=read_width),
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=f"Read ({read['flank'].capitalize()}): {temp['query_name']}<br>Read File: <b>{read_file_str}</b><br>Type: {read_type_str}<br>Range: {read['start']}-{read['end']}<br>Genomic Strand: {strand_str}"
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


def generate_known_is_locus_plot(
    chrom, is_start, is_end, is_strand, is_name,
    ref_seq, features,
    left_cov, right_cov,
    left_bam, right_bam,
    nearby_unpaired=None,
    title=""
):
    """
    Render a ±EXTENSION_PADDING window around a known IS BLASTN locus.

    Subplots (3 rows):
      1. Coverage (left+right) + GC content
      2. GenBank gene annotations
      3. Singleton flank reads (HEAD-only / TAIL-only) within the window

    The IS body is marked in every row with an amber vertical span.
    """
    if go is None:
        return None

    window_start = max(0, is_start - EXTENSION_PADDING)
    window_end   = is_end + EXTENSION_PADDING
    if ref_seq:
        window_end = min(len(ref_seq), window_end)

    is_center = (is_start + is_end) // 2

    # --- GC content ---
    gc_x, gc_y = [], []
    if ref_seq:
        gc_x, gc_y = get_gc_content_profile(ref_seq, window_start, window_end, GC_WINDOW_SIZE)

    # --- Coverage across the window ---
    l_depth = get_base_depths(left_cov,  chrom, window_start, window_end)
    r_depth = get_base_depths(right_cov, chrom, window_start, window_end)

    def _depth_to_trace(depth_map, color, name):
        coords = sorted(depth_map.keys())
        if not coords:
            return [], []
        x = [coords[0] - 1] + coords + [coords[-1] + 1]
        y = [0.0] + [depth_map[p] for p in coords] + [0.0]
        return x, y

    l_x, l_y = _depth_to_trace(l_depth, HEAD_COV_COLOR, "HEAD Coverage")
    r_x, r_y = _depth_to_trace(r_depth, TAIL_COV_COLOR, "TAIL Coverage")

    # --- Gene track ---
    overlapping_features = [f for f in (features or [])
                            if f["start"] < window_end and f["end"] > window_start]
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
    has_genes = len(overlapping_features) > 0

    # --- Singleton reads within window ---
    nearby = nearby_unpaired or []
    # Collect reads from both BAMs across the full window for singletons
    left_reads  = parse_alignments_for_flanks(left_bam,  chrom, window_start, window_end, ref_seq, "head")
    right_reads = parse_alignments_for_flanks(right_bam, chrom, window_start, window_end, ref_seq, "tail")

    # Only keep reads that belong to a nearby singleton (match by position proximity)
    nearby_positions = set()
    for h in nearby:
        try:
            pos = int(h.get('x') or h.get('y') or 0)
            nearby_positions.add(pos)
        except (TypeError, ValueError):
            pass

    # If there are nearby singletons, show all reads in the window; otherwise empty
    show_reads = left_reads + right_reads if nearby else []

    all_reads = show_reads
    templates = {}
    for r in all_reads:
        qname = r["qname"]
        if qname not in templates:
            templates[qname] = []
        templates[qname].append(r)

    template_list = []
    for qname, reads in templates.items():
        template_list.append({
            "query_name": qname,
            "reads": reads,
            "min_start": min(r["start"] for r in reads),
            "max_end": max(r["end"] for r in reads),
        })
    template_list.sort(key=lambda t: t["min_start"])
    if len(template_list) > MAX_ALIGNMENT_READS:
        step = len(template_list) / MAX_ALIGNMENT_READS
        template_list = [template_list[int(i * step)] for i in range(MAX_ALIGNMENT_READS)]

    stacked_templates = []
    row_ends = []
    for temp in template_list:
        placed = False
        for row_idx, end_coord in enumerate(row_ends):
            if end_coord + READ_PADDING < temp["min_start"]:
                row_ends[row_idx] = temp["max_end"]
                stacked_templates.append({**temp, "row": row_idx + 1})
                placed = True
                break
        if not placed:
            row_ends.append(temp["max_end"])
            stacked_templates.append({**temp, "row": len(row_ends)})
    num_read_rows = max(1, len(row_ends))

    # --- Build figure ---
    row_heights = [0.35, 0.15, 0.50] if has_genes else [0.45, 0.01, 0.54]
    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        row_heights=row_heights,
        specs=[[{"secondary_y": True}], [{"secondary_y": False}], [{"secondary_y": False}]]
    )

    # Coverage
    if l_x:
        fig.add_trace(go.Scatter(x=l_x, y=l_y, name="HEAD Coverage",
            fill='tozeroy', line=dict(color=HEAD_COV_COLOR, width=1.5),
            hovertemplate="<b>HEAD Coverage</b><br>Position: %{x}<br>Depth: %{y}<extra></extra>"),
            row=1, col=1, secondary_y=False)
    if r_x:
        fig.add_trace(go.Scatter(x=r_x, y=r_y, name="TAIL Coverage",
            fill='tozeroy', line=dict(color=TAIL_COV_COLOR, width=1.5),
            hovertemplate="<b>TAIL Coverage</b><br>Position: %{x}<br>Depth: %{y}<extra></extra>"),
            row=1, col=1, secondary_y=False)
    if gc_x:
        fig.add_trace(go.Scatter(x=gc_x, y=gc_y, name="GC Content",
            line=dict(color=GC_LINE_COLOR, width=1.2, dash='dot'),
            hovertemplate="<b>GC Content</b><br>Position: %{x}<br>GC %: %{y:.1f}%<extra></extra>"),
            row=1, col=1, secondary_y=True)

    # IS body span (amber) on all rows
    is_label = is_name.split(':')[1] if ':' in is_name else is_name
    for row_n in [1, 2, 3]:
        fig.add_vrect(
            x0=is_start, x1=is_end,
            fillcolor="rgba(234,179,8,0.12)",
            line_width=1.5, line_color="rgba(180,130,0,0.7)", line_dash="dash",
            annotation_text=is_label if row_n == 1 else "",
            annotation_position="top left",
            row=row_n, col=1
        )

    # Gene annotations
    if has_genes:
        for feat in overlapping_features:
            f_start = max(window_start, feat["start"])
            f_end   = min(window_end,   feat["end"])
            strand_str = "+" if feat["strand"] > 0 else "-"
            color = GENE_FORWARD_COLOR if feat["strand"] > 0 else GENE_REVERSE_COLOR
            label = feat["gene"] or feat["locus_tag"] or "CDS"
            hover_text = (f"Gene: {feat['gene']}<br>Locus: {feat['locus_tag']}<br>"
                          f"Product: {feat['product']}<br>Strand: {strand_str}")
            fig.add_trace(go.Scatter(
                x=[f_start, f_end], y=[feat["row"], feat["row"]],
                mode="lines+text", line=dict(color=color, width=16),
                text=[label], textposition="middle center",
                textfont=dict(color="white", size=10, weight="bold"),
                name=label, hoverinfo="text", hovertext=hover_text, showlegend=False),
                row=2, col=1)

    # Stacked reads
    for temp in stacked_templates:
        y = temp["row"]
        reads = temp["reads"]
        if len(reads) >= 2:
            r1, r2 = reads[0], reads[1]
            gap_start = min(r1["end"], r2["end"])
            gap_end   = max(r1["start"], r2["start"])
            if gap_start < gap_end:
                fig.add_trace(go.Scatter(
                    x=[gap_start, gap_end], y=[y, y], mode="lines",
                    line=dict(color="rgba(160,160,160,0.55)", width=1),
                    showlegend=False, hoverinfo="skip"), row=3, col=1)
        for read in reads:
            rc = HEAD_READ_REVERSE_COLOR if (read["flank"] == "head" and read["is_reverse"]) else \
                 HEAD_READ_FORWARD_COLOR if read["flank"] == "head" else \
                 TAIL_READ_REVERSE_COLOR if read["is_reverse"] else TAIL_READ_FORWARD_COLOR
            fig.add_trace(go.Scatter(
                x=[read["start"], read["end"]], y=[y, y], mode="lines",
                line=dict(color=rc, width=8), showlegend=False,
                hoverinfo="text",
                hovertext=f"Read ({read['flank'].capitalize()}): {temp['query_name']}<br>"
                           f"Range: {read['start']}-{read['end']}<br>"
                           f"Strand: {'Reverse' if read['is_reverse'] else 'Forward'}"),
                row=3, col=1)
            for snp in read["snps"]:
                snp_color = BASE_COLORS.get(snp["read"], "#6b7280")
                fig.add_trace(go.Scatter(
                    x=[snp["pos"]], y=[y], mode="markers+text",
                    marker=dict(color=snp_color, size=6, symbol="square"),
                    text=[snp["read"]], textposition="middle center",
                    textfont=dict(color="white", size=5, family="Courier New, monospace"),
                    showlegend=False, hoverinfo="text",
                    hovertext=f"SNP: Ref {snp['ref']} -> Read {snp['read']}<br>Position: {snp['pos']}"),
                    row=3, col=1)

    # Ref sequence (zoom-dependent)
    if ref_seq:
        ref_chars     = [ref_seq[p] for p in range(window_start, window_end)]
        ref_positions = list(range(window_start, window_end))
        fig.add_trace(go.Scatter(
            x=ref_positions, y=[0] * len(ref_positions),
            mode="text", text=ref_chars,
            textfont=dict(size=9, family="Courier New, monospace", color="#424242"),
            name="Ref Seq", showlegend=False, hoverinfo="skip",
            visible=True if (window_end - window_start <= ZOOM_THRESHOLD) else False),
            row=3, col=1)

    # Axis styling
    fig.update_xaxes(showgrid=True, gridcolor="#f3f4f6", row=1, col=1)
    fig.update_xaxes(showgrid=True, gridcolor="#f3f4f6", row=2, col=1)
    fig.update_xaxes(title_text="Genome Position", showgrid=True, gridcolor="#f3f4f6",
                     range=[window_start, window_end], row=3, col=1)
    fig.update_yaxes(title_text="Depth", row=1, col=1, secondary_y=False)
    if gc_x:
        fig.update_yaxes(title_text="GC %", range=[0, 100], row=1, col=1, secondary_y=True)
    if has_genes:
        fig.update_yaxes(showgrid=False, showticklabels=False,
                         range=[0.2, num_gene_rows + 0.8], row=2, col=1)
    else:
        fig.update_yaxes(showgrid=False, showticklabels=False, row=2, col=1)
    fig.update_yaxes(title_text="Singleton Reads", showgrid=False, showticklabels=False,
                     range=[-0.5, num_read_rows + 0.5], row=3, col=1)
    if not show_reads:
        fig.add_annotation(
            text="No singleton flanks within window",
            xref="paper", yref="paper", x=0.5, y=0.05,
            showarrow=False, font=dict(color="#9ca3af", size=13),
            row=3, col=1
        )

    fig.update_layout(
        title_text=title,
        height=PLOT_HEIGHT,
        template="plotly_white",
        margin=dict(l=60, r=60, t=50, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    return fig



def generate_report(
    table_file, report_file,
    reference_file=None,
    cutoff=6,
    known_is=None,       # list of BLASTN hits [{'chr','start','end','name','strand'}, ...]
    unpaired_hits=None,  # list of HEAD-only / TAIL-only formatted hits
    left_bam=None,       # path to left sorted BAM
    right_bam=None,      # path to right sorted BAM
):
    """
    Generates an interactive HTML report with three sections:
      1. Compact summary table (all detection classes, zero-count included)
      2. Known IS Loci (BLASTN): POSITIVE = paired flanks; NEGATIVE = singleton flanks
      3. Novel IS Loci: novel / novel (TSD) / novel (TSD)* insertions
    """
    logging.info(f"Generating HTML visualization report: {report_file} (min_y={cutoff})")

    if pd is None:
        logging.warning("Required packages for HTML report generation are missing. Skipping.")
        return

    if not os.path.exists(table_file):
        logging.warning(f"Table file {table_file} not found. Skipping report.")
        return

    df = pd.read_csv(table_file, sep='\t')

    out_dir  = os.path.dirname(os.path.abspath(table_file))
    filename = os.path.basename(table_file)

    left_bam_file  = left_bam
    right_bam_file = right_bam
    left_cov_file  = None
    right_cov_file = None

    if '__' in filename:
        sample_prefix, ref_base_with_ext = filename.split('__', 1)
        ref_base = ref_base_with_ext.replace('_table.tsv', '')
        left_cov_file  = os.path.join(out_dir, 'tmp', f"{sample_prefix}_left_{ref_base}_cov.bed")
        right_cov_file = os.path.join(out_dir, 'tmp', f"{sample_prefix}_right_{ref_base}_cov.bed")
        if not left_bam_file:
            left_bam_file  = os.path.join(out_dir, f"{sample_prefix}_left_{ref_base}.sorted.bam")
        if not right_bam_file:
            right_bam_file = os.path.join(out_dir, f"{sample_prefix}_right_{ref_base}.sorted.bam")

    left_cov  = parse_bed_cov(left_cov_file)
    right_cov = parse_bed_cov(right_cov_file)

    ref_seqs      = load_reference_seq(reference_file)
    features_dict = load_genbank_features(reference_file)

    # ================================================================
    # Helper: parse left/right position strings
    # ================================================================
    def _parse_pos(pos_str):
        try:
            a, b = map(int, str(pos_str).split('-'))
            return a, b
        except Exception:
            return None, None

    def _bam_cov_for_orient(orient):
        if orient == 'F':
            return left_bam_file, right_bam_file, left_cov, right_cov, "head", "tail"
        return right_bam_file, left_bam_file, right_cov, left_cov, "tail", "head"

    # ================================================================
    # SECTION 1 — Compact summary table
    # ================================================================
    raw_counts = {}
    if not df.empty:
        for c in df['call']:
            disp = _display_call(str(c))
            raw_counts[disp] = raw_counts.get(disp, 0) + 1
    for h in (unpaired_hits or []):
        c = _display_call(h.get('call', h.get('type', '')))
        raw_counts[c] = raw_counts.get(c, 0) + 1

    summary_rows_data = [(cls, raw_counts.get(cls, 0)) for cls in ALL_REPORT_CALL_CLASSES]
    for k, v in sorted(raw_counts.items()):
        if k not in ALL_REPORT_CALL_CLASSES:
            summary_rows_data.append((k, v))

    summary_rows_html = ""
    for cls, cnt in summary_rows_data:
        muted = ' style="color:#9ca3af;"' if cnt == 0 else ''
        summary_rows_html += f'<tr{muted}><td>{cls}</td><td style="text-align:right;font-variant-numeric:tabular-nums;">{cnt}</td></tr>\n'

    summary_html = f"""
    <table class="summary-table">
      <thead><tr><th>Detection Class</th><th style="text-align:right;">Count</th></tr></thead>
      <tbody>{summary_rows_html}</tbody>
    </table>"""

    # ================================================================
    # SECTION 2 — Known IS Loci (BLASTN)
    # ================================================================
    known_is_list = known_is or []
    known_is_section_html = ""
    known_locus_plot_ids  = []

    n_positive = 0
    n_negative = 0

    for locus_idx, is_el in enumerate(known_is_list, start=1):
        is_chrom  = is_el['chr']
        is_start  = is_el['start']
        is_end    = is_el['end']
        is_strand = is_el.get('strand', '+')
        is_name   = is_el.get('name', f'IS_{locus_idx}')
        is_center = (is_start + is_end) // 2
        is_label  = is_name.split(':')[1] if ':' in is_name else is_name

        ref_seq  = get_seq_by_chrom(ref_seqs, is_chrom)
        features = features_dict.get(is_chrom, features_dict.get(is_chrom.split('.')[0], []))

        # Find matching paired 'known'/'known*' rows within ±EXTENSION_PADDING
        matched_rows = []
        if not df.empty:
            for _, row in df.iterrows():
                if not str(row.get('call', '')).startswith('known'):
                    continue
                if row['contig'] != is_chrom:
                    continue
                try:
                    x_val = int(row['x'])
                except (TypeError, ValueError):
                    continue
                if abs(x_val - is_center) <= EXTENSION_PADDING:
                    matched_rows.append(row)

        is_positive = len(matched_rows) > 0
        if is_positive:
            n_positive += 1
            badge_class = "badge-positive"
            badge_text  = "POSITIVE"
        else:
            n_negative += 1
            badge_class = "badge-negative"
            badge_text  = "NEGATIVE"

        locus_title = (f"Locus {locus_idx}: {is_label} on {is_chrom}:{is_start}–{is_end} "
                       f"({'forward' if is_strand == '+' else 'reverse'})")

        sub_cards_html = ""

        if is_positive:
            # --- POSITIVE: one sub-card per matched paired hit ---
            for hit_idx, row in enumerate(matched_rows, start=1):
                orient = row.get('orientation', 'F')
                l_bam_src, r_bam_src, l_cov_src, r_cov_src, l_label, r_label = _bam_cov_for_orient(orient)
                l_start_h, l_end_h = _parse_pos(row['left_pos'])
                r_start_h, r_end_h = _parse_pos(row['right_pos'])
                if l_start_h is None or r_start_h is None:
                    continue
                call_disp = _display_call(row['call'])
                hit_title = (f"Hit {hit_idx}: {call_disp} ({orient}) — "
                             f"{row['left_pos']} / {row['right_pos']}")
                fig = generate_combined_alignment_plotly(
                    l_bam_src, r_bam_src, is_chrom,
                    l_start_h, l_end_h, r_start_h, r_end_h,
                    ref_seq, features,
                    l_cov_src, r_cov_src,
                    cutoff, hit_title,
                    l_label=l_label, r_label=r_label
                )
                pid = f"plot_locus{locus_idx}_hit{hit_idx}"
                known_locus_plot_ids.append(pid)
                sub_cards_html += f"""
        <div class="sub-card">
          <div class="sub-card-header">{hit_title}</div>
          <div class="card-body">{fig.to_html(full_html=False, include_plotlyjs=False, div_id=pid)}</div>
        </div>"""
        else:
            # --- NEGATIVE: singleton flanks within ±EXTENSION_PADDING ---
            nearby_unpaired = []
            for h in (unpaired_hits or []):
                if h.get('contig', '') != is_chrom:
                    continue
                try:
                    pos = int(h.get('x') or h.get('y') or 0)
                except (TypeError, ValueError):
                    continue
                if abs(pos - is_center) <= EXTENSION_PADDING:
                    nearby_unpaired.append(h)

            n_sing = len(nearby_unpaired)
            singleton_note = f" · {n_sing} singleton flank(s) nearby" if n_sing else ""
            fig = generate_known_is_locus_plot(
                chrom=is_chrom,
                is_start=is_start, is_end=is_end,
                is_strand=is_strand, is_name=is_name,
                ref_seq=ref_seq, features=features,
                left_cov=left_cov, right_cov=right_cov,
                left_bam=left_bam_file, right_bam=right_bam_file,
                nearby_unpaired=nearby_unpaired,
                title=f"No paired reads — genomic context ±{EXTENSION_PADDING:,} bp{singleton_note}",
            )
            if fig is not None:
                pid = f"plot_locus{locus_idx}_neg"
                known_locus_plot_ids.append(pid)
                sub_cards_html += f"""
        <div class="sub-card">
          <div class="sub-card-header" style="color:#92400e;">No paired flanks — singleton view{singleton_note}</div>
          <div class="card-body">{fig.to_html(full_html=False, include_plotlyjs=False, div_id=pid)}</div>
        </div>"""

        known_is_section_html += f"""
    <div class="card">
      <div class="card-header" style="display:flex;align-items:center;gap:10px;">
        {locus_title}
        <span class="badge {badge_class}">{badge_text}</span>
      </div>
      {sub_cards_html}
    </div>"""

    known_loci_intro = f"""
    <p style="color:#6b7280;margin-bottom:16px;">
      <strong>{len(known_is_list)}</strong> loci identified by BLASTN &nbsp;·&nbsp;
      <strong>{n_positive}</strong> POSITIVE (paired read evidence) &nbsp;·&nbsp;
      <strong>{n_negative}</strong> NEGATIVE (no paired reads) &nbsp;·&nbsp;
      window = ±{EXTENSION_PADDING:,} bp
    </p>""" if known_is_list else \
    "<p style='color:#9ca3af;'>No BLASTN targets file provided or no IS loci identified.</p>"

    # ================================================================
    # SECTION 3 — Novel IS Loci
    # ================================================================
    NOVEL_CLASSES   = {'novel', 'novel (TSD)', 'novel (TSD)*'}
    NOVEL_COLORS    = {
        'novel':         '#4f46e5',   # indigo
        'novel (TSD)':   '#0d9488',   # teal
        'novel (TSD)*':  '#d97706',   # amber
    }

    novel_section_html  = ""
    novel_plot_ids      = []
    novel_count         = 0

    if not df.empty:
        for idx, row in df.iterrows():
            call_raw  = str(row.get('call', ''))
            call_disp = _display_call(call_raw)
            # Strip trailing * to look up color; keep original for display
            call_base = call_raw.rstrip('*')
            if call_base not in NOVEL_CLASSES and call_raw not in NOVEL_CLASSES:
                continue

            novel_count += 1
            chrom  = row['contig']
            orient = row.get('orientation', 'F')
            accent = NOVEL_COLORS.get(call_raw, NOVEL_COLORS.get(call_base, '#4f46e5'))

            l_bam_src, r_bam_src, l_cov_src, r_cov_src, l_label, r_label = _bam_cov_for_orient(orient)
            l_start_h, l_end_h = _parse_pos(row['left_pos'])
            r_start_h, r_end_h = _parse_pos(row['right_pos'])
            if l_start_h is None or r_start_h is None:
                continue

            ref_seq  = get_seq_by_chrom(ref_seqs, chrom)
            features = features_dict.get(chrom, features_dict.get(chrom.split('.')[0], []))

            title = (f"Novel region {novel_count}: {call_disp} ({orient}) on {chrom} "
                     f"({row['left_pos']} / {row['right_pos']})")
            fig = generate_combined_alignment_plotly(
                l_bam_src, r_bam_src, chrom,
                l_start_h, l_end_h, r_start_h, r_end_h,
                ref_seq, features,
                l_cov_src, r_cov_src,
                cutoff, title,
                l_label=l_label, r_label=r_label
            )
            pid = f"plot_novel_{novel_count}"
            novel_plot_ids.append(pid)

            # Gene context note
            left_gene  = row.get('left_gene',  'N/A')
            right_gene = row.get('right_gene', 'N/A')
            gene_note  = f"{left_gene} / {right_gene}" if left_gene != 'N/A' else ""
            interrupted = row.get('gene_interruption', 'False')
            int_badge   = ' <span class="badge badge-interrupt">GENE DISRUPTED</span>' if interrupted == 'True' else ''

            novel_section_html += f"""
    <div class="card" style="border-left:4px solid {accent};">
      <div class="card-header" style="border-left:none;">
        {title}{int_badge}
        {'<span style="font-size:0.85rem;color:#6b7280;font-weight:400;margin-left:8px;">'+gene_note+'</span>' if gene_note else ''}
      </div>
      <div class="card-body">{fig.to_html(full_html=False, include_plotlyjs=False, div_id=pid)}</div>
    </div>"""

    if not novel_section_html:
        novel_section_html = "<p style='color:#9ca3af;'>No novel insertions detected.</p>"

    novel_intro = (f"<p style='color:#6b7280;margin-bottom:16px;'>"
                   f"<strong>{novel_count}</strong> novel insertion(s) — "
                   f"indigo = novel &nbsp;·&nbsp; teal = novel (TSD) &nbsp;·&nbsp; amber = novel (TSD)*</p>"
                   if novel_count else "")

    # ================================================================
    # Zoom listener JS for all plots
    # ================================================================
    all_plot_ids = known_locus_plot_ids + novel_plot_ids

    zoom_js = f"""
    <script>
    (function() {{
        const ids = {all_plot_ids};
        ids.forEach(function(pid) {{
            const gd = document.getElementById(pid);
            if (!gd) return;
            gd.on('plotly_relayout', function(ev) {{
                let x0, x1;
                ['xaxis3.range[0]','xaxis.range[0]'].forEach(function(k) {{
                    if (ev[k] !== undefined) {{ x0 = ev[k]; x1 = ev[k.replace('[0]','[1]')]; }}
                }});
                if (gd.layout && gd.layout.xaxis3 && gd.layout.xaxis3.range && x0 === undefined) {{
                    x0 = gd.layout.xaxis3.range[0]; x1 = gd.layout.xaxis3.range[1];
                }}
                if (x0 === undefined) return;
                const span = x1 - x0;
                const idxs = gd.data.map(function(d,i){{ return d.name === 'Ref Seq' ? i : -1; }}).filter(function(i){{ return i>=0; }});
                if (idxs.length) Plotly.restyle(gd, {{'visible': span <= {ZOOM_THRESHOLD}}}, idxs);
            }});
        }});
    }})();
    </script>"""

    # ================================================================
    # Full HTML assembly
    # ================================================================
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>ISLAN Mapping Report</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    body {{
      font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
      margin: 36px;
      background: #f3f4f6;
      color: #1f2937;
    }}
    .container {{
      max-width: 1400px;
      margin: auto;
      background: white;
      padding: 36px;
      border-radius: 16px;
      box-shadow: 0 8px 24px rgba(0,0,0,0.06);
    }}
    h1 {{
      font-size: 1.7rem;
      font-weight: 700;
      color: #111827;
      border-bottom: 2px solid #e5e7eb;
      padding-bottom: 14px;
      margin-top: 0;
    }}
    h2 {{
      font-size: 1.2rem;
      font-weight: 600;
      color: #374151;
      margin-top: 36px;
      margin-bottom: 10px;
    }}
    /* ── Compact summary table ─────────────────── */
    .summary-table {{
      border-collapse: collapse;
      width: auto;
      min-width: 280px;
      font-size: 0.9rem;
      margin-bottom: 28px;
    }}
    .summary-table th {{
      background: #f9fafb;
      color: #374151;
      font-weight: 600;
      padding: 8px 14px;
      border: 1px solid #e5e7eb;
      white-space: nowrap;
    }}
    .summary-table td {{
      padding: 6px 14px;
      border: 1px solid #e5e7eb;
      white-space: nowrap;
    }}
    .summary-table tr:hover td {{ background: #f0f9ff; }}
    /* ── Cards ─────────────────────────────────── */
    .card {{
      border: 1px solid #e5e7eb;
      border-radius: 12px;
      margin-bottom: 26px;
      overflow: hidden;
      box-shadow: 0 3px 8px rgba(0,0,0,0.03);
    }}
    .card-header {{
      background: #f9fafb;
      padding: 14px 20px;
      font-size: 1rem;
      font-weight: 600;
      color: #111827;
      border-bottom: 1px solid #e5e7eb;
      display: flex;
      align-items: center;
      gap: 10px;
      flex-wrap: wrap;
    }}
    .card-body {{ padding: 18px; }}
    .sub-card {{
      border-top: 1px solid #e5e7eb;
    }}
    .sub-card-header {{
      background: #fafafa;
      padding: 9px 20px;
      font-size: 0.9rem;
      font-weight: 500;
      color: #374151;
      border-bottom: 1px solid #f3f4f6;
    }}
    /* ── Badges ─────────────────────────────────── */
    .badge {{
      display: inline-block;
      padding: 3px 10px;
      border-radius: 999px;
      font-size: 0.75rem;
      font-weight: 700;
      letter-spacing: 0.05em;
      white-space: nowrap;
    }}
    .badge-positive  {{ background: #d1fae5; color: #065f46; }}
    .badge-negative  {{ background: #fef3c7; color: #92400e; }}
    .badge-interrupt {{ background: #fee2e2; color: #991b1b; }}
    /* ── FP note ─────────────────────────────────── */
    .fp-note {{
      font-size: 0.85rem;
      color: #4b5563;
      background: #f9fafb;
      padding: 10px 14px;
      border-left: 4px solid #b91c1c;
      border-radius: 6px;
      margin-bottom: 22px;
    }}
  </style>
</head>
<body>
<div class="container">

  <h1>ISLAN: Insertion Sequence Landscape Analyzer Report</h1>

  <!-- ═══ SECTION 1: Summary ═══ -->
  <h2>1 · Summary</h2>
  {summary_html}
  <div class="fp-note">
    <strong>Note on * calls:</strong> calls appended with <strong>*</strong> indicate
    cases where left and right flanking peaks overlap by more than <strong>{MAX_TSD_OVERLAP} bp</strong>
    — this is a stronger signal for a real insertion (no empty-locus gap).
  </div>

  <!-- ═══ SECTION 2: Known IS Loci ═══ -->
  <h2>2 · Known IS Loci on Reference (BLASTN)</h2>
  {known_loci_intro}
  {known_is_section_html if known_is_section_html else "<p style='color:#9ca3af;'>No known loci to display.</p>"}

  <!-- ═══ SECTION 3: Novel IS Loci ═══ -->
  <h2>3 · Novel IS Loci</h2>
  {novel_intro}
  {novel_section_html}

</div>
{zoom_js}
</body>
</html>"""

    with open(report_file, 'w') as f:
        f.write(html_content)

    logging.info("Report generated successfully.")

