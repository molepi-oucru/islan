import os
import logging
import csv

def generate_report(table_file, report_file):
    """
    Generates an interactive HTML report using Plotly and Pandas.
    """
    logging.info(f"Generating HTML visualization report: {report_file}")
    try:
        import pandas as pd
        import plotly.express as px
        import plotly.graph_objects as go
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

    data = []
    for _, row in df.iterrows():
        hit_type = row.get('call', 'unknown')
        x_val = row.get('x', 'N/A')
        y_val = row.get('y', 'N/A')
        gap_val = row.get('gap', -1)
        
        l_depth = float(row['left_depth_median']) if row.get('left_depth_median', 'N/A') != 'N/A' else 0
        r_depth = float(row['right_depth_median']) if row.get('right_depth_median', 'N/A') != 'N/A' else 0
        
        try:
            gap = int(gap_val) if pd.notna(gap_val) and gap_val != 'N/A' else -1
        except ValueError:
            gap = -1
            
        data.append({
            'type': hit_type,
            'max_depth': max(l_depth, r_depth),
            'min_depth': min(l_depth, r_depth) if min(l_depth, r_depth) > 0 else max(l_depth, r_depth),
            'gap': gap,
            'left_pos': str(x_val),
            'right_pos': str(y_val),
            'gene': row.get('left_gene', '')
        })

    pdf = pd.DataFrame(data)

    # 1. Summary Table
    summary = pdf.groupby('type').size().reset_index(name='count')
    summary_html = summary.to_html(classes='table table-striped', index=False)

    # 2. Scatter Plot: Gap vs Depth (Only for Pairs)
    pairs = pdf[pdf['gap'] >= 0]
    scatter_html = ""
    if not pairs.empty:
        fig1 = px.scatter(
            pairs, x='gap', y='min_depth', color='type',
            hover_data=['left_pos', 'right_pos', 'gene'],
            log_y=True,
            title='Gap Distance vs Minimum Pair Depth',
            labels={'gap': 'Gap Distance (bp)', 'min_depth': 'Minimum Depth of Pair (Log Scale)'}
        )
        fig1.add_vline(x=100, line_dash="dash", line_color="red", annotation_text="100bp Gap Limit")
        scatter_html = fig1.to_html(full_html=False, include_plotlyjs='cdn')

    # 3. Histogram of Depths
    fig2 = px.histogram(
        pdf, x='max_depth', color='type',
        log_x=True,
        title='Distribution of Maximum Peak Depths',
        labels={'max_depth': 'Max Depth (Log Scale)'}
    )
    hist_html = fig2.to_html(full_html=False, include_plotlyjs='cdn')

    html_content = f"""
    <html>
    <head>
        <title>IS-Seq Mapping Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .container {{ max-width: 1200px; margin: auto; }}
            .table {{ border-collapse: collapse; width: 50%; margin-bottom: 20px; }}
            .table th, .table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .table th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>IS-Seq Mapping Report</h1>
            <h2>Summary Statistics</h2>
            {summary_html}
            <h2>Visualization</h2>
            {scatter_html}
            <br>
            {hist_html}
        </div>
    </body>
    </html>
    """

    with open(report_file, 'w') as f:
        f.write(html_content)
        
    logging.info(f"Report generated successfully.")
