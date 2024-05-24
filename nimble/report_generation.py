import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Environment, FileSystemLoader

def classify_max_scores(row):
    filter_threshold = 67

    # Replace 0 scores with NaN for max calculation
    scores = {'r1_forward_score': row['r1_forward_score'],
              'r1_reverse_score': row['r1_reverse_score'],
              'r2_forward_score': row['r2_forward_score'],
              'r2_reverse_score': row['r2_reverse_score']}
    for key in scores:
        if scores[key] == 0:
            scores[key] = np.nan
    
    # Calculate max scores for r1 and r2
    max_r1 = max(scores['r1_forward_score'], scores['r1_reverse_score'])
    max_r2 = max(scores['r2_forward_score'], scores['r2_reverse_score'])
    
    # Apply filter conditions
    filtered = (np.isnan(max_r1) or max_r1 >= filter_threshold) and (np.isnan(max_r2) or max_r2 >= filter_threshold)
    return 'unfiltered' if filtered else 'filtered'

def classify_scores_with_filter(row):
    threshold = 50
    scores = [row['r1_forward_score'], row['r1_reverse_score'], row['r2_forward_score'], row['r2_reverse_score']]
    filters = [row['r1_filter_forward'], row['r1_filter_reverse'], row['r2_filter_forward'], row['r2_filter_reverse']]
    score_columns = ['r1_forward_score', 'r1_reverse_score', 'r2_forward_score', 'r2_reverse_score']

    # Check if any scores are above the threshold
    score_names = [score_columns[i] for i, score in enumerate(scores) if score > threshold]
    score_label = '+'.join(score_names) if score_names else 'None'

    # Check filters
    filter_label = 'score_below_threshold' if any(f == 'Score Below Threshold' for f in filters) else 'no_score_below_threshold'
    
    # Combine labels
    combined_label = score_label + '+' + filter_label
    return combined_label

def tag_dataframe(df):
    df['Max_Score_Filter'] = df.apply(classify_max_scores, axis=1)
    df['Classification'] = df.apply(classify_scores_with_filter, axis=1)
    df['Is_Filtered'] = df.apply(lambda row: 'unfiltered' if row['Classification'] in [
        'None+no_score_below_threshold', 'None+score_below_threshold', 
        'r1_reverse_score+no_score_below_threshold', 'r1_forward_score+no_score_below_threshold',
        'r2_forward_score+no_score_below_threshold', 'r2_reverse_score+no_score_below_threshold'
    ] or row['Max_Score_Filter'] == 'filtered' else 'filtered', axis=1)
    return df

def filter_dataframe(df):
    df = tag_dataframe(df)
    return df[df['Is_Filtered'] == 'unfiltered']


def generate_plots_for_feature(df, output_dir, nimble_feature):
    # Number of unique values
    num_UMIs = df['r1_UB'].nunique()
    num_cells = df['r1_CB'].nunique()
    
    # Density graph for r1_POS and r2_POS
    plt.figure(figsize=(12, 6))
    sns.kdeplot(df[df['r1_POS'] > 0]['r1_POS'], bw_adjust=0.5, label='r1_POS', color='blue', linestyle='--')
    sns.kdeplot(df[df['r2_POS'] > 0]['r2_POS'], bw_adjust=0.5, label='r2_POS', color='red', linestyle='-')
    plt.title('Density Plot of Positions')
    plt.legend()
    plt.xlim(df.loc[df['r1_POS'] > 0, 'r1_POS'].min(), df.loc[df['r2_POS'] > 0, 'r2_POS'].max())
    plt.xlabel('Position')
    plt.savefig(f"{output_dir}/{nimble_feature}_positions_density.png")
    plt.close()

    # Violin plot for scores
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=df[['r1_forward_score', 'r1_reverse_score', 'r2_forward_score', 'r2_reverse_score']])
    plt.title('Violin Plot of Scores')
    plt.savefig(f"{output_dir}/{nimble_feature}_scores_violin.png")
    plt.close()

    # Bar plot for Classifications
    plt.figure(figsize=(12, 6))
    ax = df['Classification'].value_counts().plot(kind='bar')
    plt.title('Bar Plot of Classifications')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{nimble_feature}_classifications_bar.png")
    plt.close()

    # Load Jinja2 template
    current_dir = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(current_dir, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template('feature_report_template.html')
    
    # Render report
    html_content = template.render(num_UMIs=num_UMIs, num_cells=num_cells, nimble_feature=nimble_feature,
                                   pos_density_image=f"{nimble_feature}_positions_density.png",
                                   scores_violin_image=f"{nimble_feature}_scores_violin.png",
                                   classifications_bar_image=f"{nimble_feature}_classifications_bar.png")
    with open(f"{output_dir}/{nimble_feature}_report.html", 'w') as f:
        f.write(html_content)

def generate_plots(df, output_dir):
    df = tag_dataframe(df)
    valid_features = [feature for feature in df['nimble_features'].dropna().unique() if ',' not in feature]
    report_files = []

    for feature in valid_features:
        print(f"Writing plots for feature {feature}")
        feature_df = df[df['nimble_features'] == feature]
        report_filename = f"{feature}_report.html"
        report_file_path = os.path.join(output_dir, report_filename)
        generate_plots_for_feature(feature_df, output_dir, feature)
        report_files.append(report_file_path)

    print("Writing final report")
    concatenate_reports(report_files, output_dir)

def concatenate_reports(report_files, output_dir):
    main_report_path = os.path.join(output_dir, "final_report.html")
    with open(main_report_path, 'w') as main_file:
        main_file.write("""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Combined Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                h1, h2 { color: #333; }
                img { width: 100%; height: auto; margin-top: 20px; }
                .statistics { margin-top: 20px; }
                .statistics p { font-size: 16px; color: #666; }
            </style>
        </head>
        <body>
        """)

        for report_file in report_files:
            with open(report_file, 'r') as file:
                main_file.write(file.read() + "\n<br/>\n")
        
        main_file.write("</body></html>")