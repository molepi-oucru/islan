import os
import logging
import csv

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

def get_depth_profile(cov, chrom, start, end, cutoff):
    coords = list(range(start, end))
    depth_map = {c: 0.0 for c in coords}
    
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
            if pos in depth_map:
                depth_map[pos] = c_depth
                
    sorted_coords = sorted(depth_map.keys())
    sorted_depths = [max(float(cutoff), depth_map[c]) for c in sorted_coords]
    return sorted_coords, sorted_depths

def generate_report(table_file, report_file, cutoff=6):
    """
    Generates an interactive HTML report using Plotly and Pandas.
    """
    logging.info(f"Generating HTML visualization report: {report_file} (min_y={cutoff})")
    try:
        import pandas as pd
        import numpy as np
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        logging.warning("Pandas or Plotly not installed. Skipping HTML report generation. Please run: pip install pandas plotly")
        return

    if not os.path.exists(table_file):
        logging.warning(f"Table file {table_file} not found. Skipping report.")
        return

    df = pd.read_csv(table_file, sep='\t')
    
    if df.empty:
        logging.warning("Table is empty. Skipping report.")
        return

    # Reconstruct bedGraph paths from table_file path
    out_dir = os.path.dirname(os.path.abspath(table_file))
    filename = os.path.basename(table_file)
    
    left_cov_file = None
    right_cov_file = None
    if '__' in filename:
        sample_prefix, ref_base_with_ext = filename.split('__', 1)
        ref_base = ref_base_with_ext.replace('_table.tsv', '')
        left_cov_file = os.path.join(out_dir, 'tmp', f"{sample_prefix}_left_{ref_base}_cov.bed")
        right_cov_file = os.path.join(out_dir, 'tmp', f"{sample_prefix}_right_{ref_base}_cov.bed")

    # Load coverages
    left_cov = parse_bed_cov(left_cov_file)
    right_cov = parse_bed_cov(right_cov_file)

    # 1. Summary Table HTML
    summary = df.groupby('call').size().reset_index(name='count')
    summary_html = summary.to_html(classes='table table-striped', index=False)

    # 2. Generate overlaid flanking coverage plots for each region
    plots_html = ""
    if not left_cov and not right_cov:
        plots_html = "<p class='warning'>Raw coverage files not found. Flanking coverage plots could not be generated.</p>"
    else:
        # Create subplots grid: N rows, 1 col
        num_regions = len(df)
        fig = make_subplots(
            rows=num_regions, cols=1,
            subplot_titles=[f"Region {idx+1} ({row['call']}) Flanking Coverage" for idx, row in df.iterrows()],
            vertical_spacing=max(0.02, 0.2 / num_regions) if num_regions > 1 else 0.05
        )
        
        y_range = [np.log10(float(cutoff)), np.log10(30000.0)]
        show_head_legend = True
        show_tail_legend = True

        for idx, row in df.iterrows():
            row_num = idx + 1
            chrom = row['contig']
            orient = row.get('orientation', 'F')
            region_name = row.get('region', f"region_{row_num}")
            
            try:
                l_start, l_end = map(int, str(row['left_pos']).split('-'))
                r_start, r_end = map(int, str(row['right_pos']).split('-'))
            except Exception as e:
                logging.warning(f"Could not parse flank coordinates for row {row_num}: {e}")
                continue

            if orient == 'F':
                l_coords, l_depths = get_depth_profile(left_cov, chrom, l_start, l_end, cutoff)
                r_coords, r_depths = get_depth_profile(right_cov, chrom, r_start, r_end, cutoff)
                l_is_head = True
                l_show = show_head_legend
                if show_head_legend: show_head_legend = False
                
                r_is_head = False
                r_show = show_tail_legend
                if show_tail_legend: show_tail_legend = False
            else:
                l_coords, l_depths = get_depth_profile(right_cov, chrom, l_start, l_end, cutoff)
                r_coords, r_depths = get_depth_profile(left_cov, chrom, r_start, r_end, cutoff)
                l_is_head = False
                l_show = show_tail_legend
                if show_tail_legend: show_tail_legend = False
                
                r_is_head = True
                r_show = show_head_legend
                if show_head_legend: show_head_legend = False

            # Add Left Flank Trace
            fig.add_trace(
                go.Scatter(
                    x=l_coords, y=l_depths,
                    name="HEAD" if l_is_head else "TAIL",
                    legendgroup="HEAD" if l_is_head else "TAIL",
                    showlegend=l_show,
                    mode='lines+markers',
                    line=dict(color='royalblue' if l_is_head else 'firebrick', width=1.5),
                    marker=dict(size=4),
                    hovertemplate=f"<b>{region_name} {'HEAD' if l_is_head else 'TAIL'}</b><br>Pos: %{{x}}<br>Depth: %{{y}}<extra></extra>"
                ),
                row=row_num, col=1
            )
            
            # Add Right Flank Trace
            fig.add_trace(
                go.Scatter(
                    x=r_coords, y=r_depths,
                    name="HEAD" if r_is_head else "TAIL",
                    legendgroup="HEAD" if r_is_head else "TAIL",
                    showlegend=r_show,
                    mode='lines+markers',
                    line=dict(color='royalblue' if r_is_head else 'firebrick', width=1.5),
                    marker=dict(size=4),
                    hovertemplate=f"<b>{region_name} {'HEAD' if r_is_head else 'TAIL'}</b><br>Pos: %{{x}}<br>Depth: %{{y}}<extra></extra>"
                ),
                row=row_num, col=1
            )
            
            # Formatting for this subplot
            fig.update_xaxes(title_text="Genome Position", row=row_num, col=1)
            fig.update_yaxes(
                title_text="Depth (Log Scale)",
                type="log",
                range=y_range,
                dtick=1.0,
                row=row_num, col=1
            )

        # Dynamic figure height based on number of subplots
        fig_height = max(400, 350 * num_regions)
        fig.update_layout(
            title_text=f"Overlaid Flanking Coverage Profiles (Log-Scale Depth, min_y={cutoff})",
            height=fig_height,
            width=1000,
            showlegend=True,
            template="plotly_white"
        )
        
        plots_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IS-Seq Mapping Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 40px;
                background-color: #f8f9fa;
                color: #333;
            }}
            .container {{
                max-width: 1100px;
                margin: auto;
                background: white;
                padding: 30px;
                border-radius: 12px;
                box-shadow: 0 4px 15px rgba(0,0,0,0.05);
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 2px solid #ecf0f1;
                padding-bottom: 15px;
                font-weight: 600;
            }}
            h2 {{
                color: #34495e;
                margin-top: 30px;
                font-weight: 500;
            }}
            .table {{
                border-collapse: collapse;
                width: 100%;
                margin-top: 15px;
                margin-bottom: 30px;
            }}
            .table th, .table td {{
                border: 1px solid #e2e8f0;
                padding: 12px 15px;
                text-align: left;
            }}
            .table th {{
                background-color: #f7fafc;
                color: #4a5568;
                font-weight: 600;
            }}
            .table tr:nth-child(even) {{
                background-color: #fcfcfc;
            }}
            .warning {{
                color: #e53e3e;
                background: #fff5f5;
                padding: 12px;
                border-left: 4px solid #e53e3e;
                border-radius: 4px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>IS-Seq Mapping Report</h1>
            <h2>Summary of Detected Insertions</h2>
            {summary_html}
            <h2>Flanking Region Coverage Depth Profiles</h2>
            {plots_html}
        </div>
    </body>
    </html>
    """

    with open(report_file, 'w') as f:
        f.write(html_content)
        
    logging.info(f"Report generated successfully.")
