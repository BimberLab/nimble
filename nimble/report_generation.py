import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Environment, FileSystemLoader
import base64
from io import BytesIO

def position_density(df, threshold=150000, boundary_expansion=10000):
    # Combine, sort, and get unique positions
    all_positions = np.sort(np.unique(np.concatenate((df['r1_POS'].dropna().values, df['r2_POS'].dropna().values))))

    # Detect discontinuities
    diffs = np.diff(all_positions)
    break_indices = np.where(diffs > threshold)[0]

    # Adjust threshold if there are more than 5 plots
    while len(break_indices) > 4:
        threshold *= 2
        break_indices = np.where(diffs > threshold)[0]

    if len(break_indices) == 0:
        # If no discontinuities, create a normal axis
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.kdeplot(df['r1_POS'], bw_adjust=0.5, label='r1_POS', color='blue', linestyle='--', ax=ax, warn_singular=False)
        sns.kdeplot(df['r2_POS'], bw_adjust=0.5, label='r2_POS', color='red', linestyle='-', ax=ax, warn_singular=False)
        ax.legend()
        ax.set_xlim(min(df['r1_POS'].min(), df['r2_POS'].min()), max(df['r1_POS'].max(), df['r2_POS'].max()))
        ax.set_xlabel('Position within genome')
        ax.set_ylabel('Density')
        pos_density_image = generate_base64_image(fig)
        plt.close(fig)
        return pos_density_image

    # Create multiple segments for each region between breaks
    segments = []
    prev_index = 0
    for idx in break_indices:
        start = max(0, all_positions[prev_index] - boundary_expansion)
        end = all_positions[idx] + boundary_expansion
        segments.append((start, end))
        prev_index = idx + 1
    start = max(0, all_positions[break_indices[-1] + 1] - boundary_expansion)
    end = all_positions.max() + boundary_expansion
    segments.append((start, end))

    # Filter out empty segments
    valid_segments = []
    for start, end in segments:
        if len(df['r1_POS'][(df['r1_POS'] >= start) & (df['r1_POS'] <= end)]) > 1 or \
           len(df['r2_POS'][(df['r2_POS'] >= start) & (df['r2_POS'] <= end)]) > 1:
            valid_segments.append((start, end))

    # Plot with discontinuities using subplots
    num_segments = len(valid_segments)
    fig, axes = plt.subplots(1, num_segments, figsize=(14, 6), sharey=True)

    if num_segments == 1:
        axes = [axes]  # To make it iterable if there's only one segment

    for i, (start, end) in enumerate(valid_segments):
        sns.kdeplot(df['r1_POS'][(df['r1_POS'] >= start) & (df['r1_POS'] <= end)], 
                    bw_adjust=0.5, color='blue', linestyle='--', ax=axes[i], warn_singular=False)
        sns.kdeplot(df['r2_POS'][(df['r2_POS'] >= start) & (df['r2_POS'] <= end)], 
                    bw_adjust=0.5, color='red', linestyle='-', ax=axes[i], warn_singular=False)
        axes[i].set_xlim(start, end)
        axes[i].tick_params(left=True)
        axes[i].tick_params(bottom=True)
        axes[i].set_xlabel('')
        axes[i].set_ylabel('')

    fig.text(0.5, 0.01, 'Position within genome', ha='center')
    fig.text(0.01, 0.5, 'Density', va='center', rotation='vertical')

    handles = [plt.Line2D([0,1],[0,1], color='blue', linestyle='--', label='r1_POS'),
               plt.Line2D([0,1],[0,1], color='red', linestyle='-', label='r2_POS')]
    axes[-1].legend(handles=handles, loc='upper right')

    pos_density_image = generate_base64_image(fig)
    plt.close(fig)
    return pos_density_image

# Violin plot for scores
def score_violin(df):
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(data=df[['r1_forward_score', 'r2_forward_score']], ax=ax)
    
    # Set fixed tick positions and labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['r1 score in bp', 'r2 score in bp'])
    
    scores_violin_image = generate_base64_image(fig)
    plt.close(fig)
    return scores_violin_image

def feature_confusion_matrix(df):
    df.loc[df['r1_GN'].isnull(), 'r1_GN'] = 'Not called'
    
    feature_counts = df.groupby(['nimble_features', 'r1_GN']).size().unstack(fill_value=0)
    
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.heatmap(feature_counts, annot=True, fmt='d', cmap='viridis', ax=ax)
    ax.set_xlabel('Input file call')
    ax.set_ylabel('Nimble call')
    feature_conf_matrix_image = generate_base64_image(fig)
    plt.close(fig)
    return feature_conf_matrix_image


def generate_plots_for_feature(df, nimble_feature):
    # Calculate coarse statistics that get reported at the top
    num_UMIs = df['r1_UB'].nunique()
    num_cells = df['r1_CB'].nunique()
    
    pos_density_image = position_density(df)
    scores_violin_image = score_violin(df)
    feature_conf_matrix_image = feature_confusion_matrix(df)

    # Load Jinja2 template
    current_dir = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(current_dir, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template('feature_report_template.html')
    
    # Render report
    html_content = template.render(num_UMIs=num_UMIs, num_cells=num_cells, nimble_feature=nimble_feature,
                                   pos_density_image=pos_density_image,
                                   scores_violin_image=scores_violin_image,
                                   feature_conf_matrix_image=feature_conf_matrix_image)
    return html_content

def generate_plots(df, output_file):
    valid_features = [feature for feature in df['nimble_features'].dropna().unique() if ',' not in feature]
    reports = []

    for feature in valid_features:
        print(f"Generating plots for feature {feature}")
        feature_df = df[df['nimble_features'] == feature]
        report_content = generate_plots_for_feature(feature_df, feature)
        reports.append(report_content)

    print("Writing final report")
    concatenate_reports(reports, output_file)

def generate_base64_image(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    return img_base64

def concatenate_reports(reports, output_file):
    with open(output_file, 'w') as f:
        f.write("""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Combined Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                h1, h2 { color: #333; }
                img { width: 70%; height: 50%; margin-top: 20px; }
            </style>
        </head>
        <body>
        """)

        for report in reports:
            f.write(report + "\n<br/>\n")
        
        f.write("</body></html>")
