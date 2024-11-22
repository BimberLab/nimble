import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Environment, FileSystemLoader
import base64
from io import BytesIO
from nimble.utils import per_umi_thresholding, umi_intersection

def plot_top_two_feature_ratio(df):
    # Calculate feature scores per UMI
    def get_top_two_ratio(umi_group):
        counts = []
        total_score = 0
        for _, row in umi_group.iterrows():
            features = row['features'].split(',')
            score_per_feature = row['nimble_score'] / len(features)
            total_score += row['nimble_score']
            for feature in features:
                counts.append({'feature': feature, 'score': score_per_feature})
        feature_scores = pd.DataFrame(counts).groupby('feature')['score'].sum()
        top_two_scores = feature_scores.nlargest(2).values
        if len(top_two_scores) == 2:
            ratio = top_two_scores[0] / (top_two_scores[0] + top_two_scores[1])
        elif len(top_two_scores) == 1:
            ratio = 1.0  # Only one feature present
        else:
            ratio = np.nan  # No features present
        return ratio

    # Apply to each UMI
    df['top_two_ratio'] = df.groupby(['cb', 'umi']).apply(get_top_two_ratio).values

    # Remove NaN values
    df_valid = df.dropna(subset=['top_two_ratio'])

    # Plot density plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.kdeplot(df_valid['top_two_ratio'], ax=ax)
    ax.set_title('Density Plot of Top Feature Ratio over Top Two Features per UMI')
    ax.set_xlabel('Top Feature Ratio')
    ax.set_ylabel('Density')
    plt.show()

def plot_feature_score_distributions(df, top_n=20):
    # Expand features and scores per row
    def explode_features(row):
        features = row['features'].split(',')
        score_per_feature = row['nimble_score'] / len(features)
        return pd.DataFrame({
            'cb': row['cb'],
            'umi': row['umi'],
            'feature': features,
            'score': [score_per_feature]*len(features)
        })

    df_exploded = pd.concat(df.apply(explode_features, axis=1).tolist(), ignore_index=True)

    # Identify top N features
    top_features = df_exploded['feature'].value_counts().nlargest(top_n).index.tolist()

    # Filter for top features
    df_top_features = df_exploded[df_exploded['feature'].isin(top_features)]

    # Plot violin plots
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(data=df_top_features, x='feature', y='score', ax=ax, order=top_features)
    ax.set_title(f'Score Distributions for Top {top_n} Features')
    ax.set_xlabel('Feature')
    ax.set_ylabel('Score')
    plt.xticks(rotation=45)
    plt.show()

def plot_stacked_bar_feature_ratios(df, sample_size=50):
    # Sample UMIs
    sampled_umis = df[['cb', 'umi']].drop_duplicates().sample(n=sample_size, random_state=42)

    # Filter data for sampled UMIs
    df_sampled = pd.merge(df, sampled_umis, on=['cb', 'umi'])

    # Calculate feature ratios per UMI
    def calculate_feature_ratios(umi_group):
        counts = []
        for _, row in umi_group.iterrows():
            features = row['features'].split(',')
            score_per_feature = row['nimble_score'] / len(features)
            for feature in features:
                counts.append({'feature': feature, 'score': score_per_feature})
        feature_scores = pd.DataFrame(counts).groupby('feature')['score'].sum()
        total_score = feature_scores.sum()
        feature_ratios = feature_scores / total_score
        return feature_ratios.reset_index()

    umi_feature_ratios = df_sampled.groupby(['cb', 'umi']).apply(calculate_feature_ratios).reset_index()
    umi_feature_ratios.rename(columns={'level_2': 'index'}, inplace=True)

    # Pivot data for plotting
    plot_data = umi_feature_ratios.pivot_table(
        index=['cb', 'umi'],
        columns='feature',
        values='score',
        fill_value=0
    )

    # Plot stacked bar chart
    plot_data.plot(kind='bar', stacked=True, figsize=(15, 8), legend=False)
    plt.title('Feature Ratios per UMI (Sampled)')
    plt.xlabel('UMI (cb, umi)')
    plt.ylabel('Feature Ratio')
    plt.legend(title='Feature', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

def plot_feature_heatmap(df, top_n_features=20, sample_size=100):
    # Expand features and scores per row
    def explode_features(row):
        features = row['features'].split(',')
        score_per_feature = row['nimble_score'] / len(features)
        return pd.DataFrame({
            'umi': row['umi'],
            'feature': features,
            'score': [score_per_feature]*len(features)
        })

    df_exploded = pd.concat(df.apply(explode_features, axis=1).tolist(), ignore_index=True)

    # Identify top N features
    top_features = df_exploded['feature'].value_counts().nlargest(top_n_features).index.tolist()

    # Filter for top features
    df_top_features = df_exploded[df_exploded['feature'].isin(top_features)]

    # Sample UMIs
    sampled_umis = df_top_features['umi'].drop_duplicates().sample(n=sample_size, random_state=42)

    # Filter data for sampled UMIs
    df_sampled = df_top_features[df_top_features['umi'].isin(sampled_umis)]

    # Pivot data to create matrix for heatmap
    heatmap_data = df_sampled.pivot_table(
        index='umi',
        columns='feature',
        values='score',
        aggfunc='sum',
        fill_value=0
    )

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_data, cmap='viridis')
    plt.title('Feature Scores Heatmap (Sampled UMIs and Top Features)')
    plt.xlabel('Feature')
    plt.ylabel('UMI')
    plt.show()

    # Expand features and scores per row
    def explode_features(row):
        features = row['features'].split(',')
        score_per_feature = row['nimble_score'] / len(features)
        return pd.DataFrame({
            'cb': row['cb'],
            'umi': row['umi'],
            'feature': features,
            'score': [score_per_feature]*len(features)
        })

    df_exploded = pd.concat(df.apply(explode_features, axis=1).tolist(), ignore_index=True)

    # Calculate total score per UMI
    umi_total_scores = df_exploded.groupby(['cb', 'umi'])['score'].sum().reset_index(name='total_score')

    # Merge total scores back to exploded data
    df_exploded = pd.merge(df_exploded, umi_total_scores, on=['cb', 'umi'])

    # Calculate feature ratios per UMI
    df_exploded['feature_ratio'] = df_exploded['score'] / df_exploded['total_score']

    # Identify top N features
    top_features = df_exploded['feature'].value_counts().nlargest(top_n).index.tolist()

    # Filter for top features
    df_top_features = df_exploded[df_exploded['feature'].isin(top_features)]

    # Plot boxplots
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(data=df_top_features, x='feature', y='feature_ratio', ax=ax, order=top_features)
    ax.set_title(f'Top Feature Ratios Across UMIs for Top {top_n} Features')
    ax.set_xlabel('Feature')
    ax.set_ylabel('Feature Ratio')
    plt.xticks(rotation=45)
    plt.show()

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
    
    if num_segments == 0:
        raise ValueError("No valid segments found for position density plot.")

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
    try:
        # Calculate coarse statistics that get reported at the top
        num_UMIs = df['r1_UB'].nunique()
        num_cells = df['r1_CB'].nunique()

        # Compute the number of read-mates expressing the feature per UMI
        # First, group by UMI and count the number of read-mates
        umi_read_mate_counts = df.groupby('r1_UB').size().reset_index(name='read_mate_count')

        # Plot the distribution of read-mate counts per UMI
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.histplot(umi_read_mate_counts['read_mate_count'], bins=range(1, umi_read_mate_counts['read_mate_count'].max() + 2), ax=ax, discrete=True)
        ax.set_title(f'Distribution of Read-Mates per UMI for Feature {nimble_feature}')
        ax.set_xlabel('Number of Read-Mates per UMI')
        ax.set_ylabel('Frequency')
        read_mates_per_umi_image = generate_base64_image(fig)
        plt.close(fig)
        
        pos_density_image = position_density(df)
        scores_violin_image = score_violin(df)
        feature_conf_matrix_image = feature_confusion_matrix(df)

        # Load Jinja2 template
        current_dir = os.path.dirname(os.path.abspath(__file__))
        templates_dir = os.path.join(current_dir, 'templates')
        env = Environment(loader=FileSystemLoader(templates_dir))
        template = env.get_template('feature_report_template.html')
        
        # Render report
        html_content = template.render(
            num_UMIs=num_UMIs,
            num_cells=num_cells,
            nimble_feature=nimble_feature,
            read_mates_per_umi_image=read_mates_per_umi_image,
            pos_density_image=pos_density_image,
            scores_violin_image=scores_violin_image,
            feature_conf_matrix_image=feature_conf_matrix_image
        )
        return html_content
    except Exception as e:
        print(f"Error generating plots for feature {nimble_feature}: {e}")
        return f"<h2>Error generating plots for feature {nimble_feature}</h2><p>{e}</p>"

def generate_plots(df, output_file):
    reports = []

    # First, generate the per-UMI intersection pages
    umi_reports = generate_umi_intersection_reports(df)
    reports.extend(umi_reports)

    # Then, proceed with the per-feature plots as before
    valid_features = [feature for feature in df['nimble_features'].dropna().unique() if ',' not in feature]

    for feature in valid_features:
        print(f"Generating plots for feature {feature}")
        feature_df = df[df['nimble_features'] == feature]
        report_content = generate_plots_for_feature(feature_df, feature)
        reports.append(report_content)

    print("Writing final report")
    concatenate_reports(reports, output_file)

def generate_umi_intersection_reports(df):
    print("Computing summary plots for whole dataset")
    # Similar preprocessing to __main__ report() so UMI intersection works
    df = df.rename(columns={
        'r1_CB': 'cb',
        'r1_UB': 'umi',
        'nimble_features': 'features'
    })

    df = df.dropna(subset=['features', 'umi', 'cb', 'nimble_score'])
    df = df[(df['features'] != '') & (df['umi'] != '') & (df['cb'] != '')]

    # Ensure features are sorted
    df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))
    df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()

    # Generate summary plots before thresholding and intersection
    # Number of UMIs total
    total_umis = df[['cb', 'umi']].drop_duplicates().shape[0]

    # Number of cells total
    total_cells = df['cb'].nunique()

    # Distribution of the number of read-mates per UMI
    umi_read_counts = df.groupby(['cb', 'umi']).size().reset_index(name='read_mate_count')
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(umi_read_counts['read_mate_count'], bins=50, ax=ax)
    ax.set_title('Distribution of Number of Read-Mates per UMI')
    ax.set_xlabel('Number of Read-Mates per UMI')
    ax.set_ylabel('Frequency')
    umi_read_counts_image = generate_base64_image(fig)
    plt.close(fig)

    # Distribution of the number of features per UMI (integer counts)
    # Calculate the number of unique features per UMI
    umi_feature_counts = df.groupby(['cb', 'umi'])['features'].apply(
        lambda x: len(set(','.join(x).split(',')))
    ).reset_index(name='num_features')

    fig, ax = plt.subplots(figsize=(8, 6))
    max_num_features = int(umi_feature_counts['num_features'].max())
    sns.histplot(
        umi_feature_counts['num_features'],
        bins=range(1, max_num_features + 2),
        ax=ax,
        discrete=True
    )
    ax.set_title('Distribution of Number of Features per UMI')
    ax.set_xlabel('Number of Features per UMI')
    ax.set_ylabel('Frequency')
    umi_feature_counts_image = generate_base64_image(fig)
    plt.close(fig)

    # Distribution of the number of read-mates per cell
    cell_read_counts = df.groupby(['cb']).size().reset_index(name='read_mate_count')
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(cell_read_counts['read_mate_count'], bins=50, ax=ax)
    ax.set_title('Distribution of Number of Read-Mates per Cell')
    ax.set_xlabel('Number of Read-Mates per Cell')
    ax.set_ylabel('Frequency')
    cell_read_counts_image = generate_base64_image(fig)
    plt.close(fig)

    # Distribution of the number of unique features per cell (integer counts)
    cell_feature_counts = df.groupby('cb')['features'].apply(
        lambda x: len(set(','.join(x).split(',')))
    ).reset_index(name='num_features')

    fig, ax = plt.subplots(figsize=(8, 6))
    max_num_features = int(cell_feature_counts['num_features'].max())
    sns.histplot(
        cell_feature_counts['num_features'],
        bins=range(1, max_num_features + 2),
        ax=ax,
        discrete=True
    )
    ax.set_title('Distribution of Number of Unique Features per Cell')
    ax.set_xlabel('Number of Features per Cell')
    ax.set_ylabel('Frequency')
    cell_feature_counts_image = generate_base64_image(fig)
    plt.close(fig)

    # Additional plots for the ratio of top two features per UMI
    def get_top_two_ratio(umi_group):
        counts = []
        total_score = 0
        for _, row in umi_group.iterrows():
            features = row['features'].split(',')
            score_per_feature = row['nimble_score'] / len(features)
            total_score += row['nimble_score']
            for feature in features:
                counts.append({'feature': feature, 'score': score_per_feature})
        feature_scores = pd.DataFrame(counts).groupby('feature')['score'].sum()
        top_two_scores = feature_scores.nlargest(2).values
        if len(top_two_scores) == 2:
            ratio = top_two_scores[0] / (top_two_scores[0] + top_two_scores[1])
        elif len(top_two_scores) == 1:
            ratio = 1.0  # Only one feature present
        else:
            ratio = np.nan  # No features present
        return ratio

    umi_read_counts['top_two_ratio'] = df.groupby(['cb', 'umi']).apply(get_top_two_ratio).values

    umi_ratio_valid = umi_read_counts.dropna(subset=['top_two_ratio'])

    # Plot density plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.kdeplot(umi_ratio_valid['top_two_ratio'], ax=ax)
    ax.set_title('Density Plot of Top Feature Ratio over Top Two Features per UMI')
    ax.set_xlabel('Top Feature Ratio')
    ax.set_ylabel('Density')
    top_two_ratio_image = generate_base64_image(fig)
    plt.close(fig)

    # Prepare variables for the template
    summary_stats = {
        'total_umis': total_umis,
        'total_cells': total_cells
    }

    # threshold simulations
    thresholds = [0, 0.05, 0.10, 0.20, 0.50]  # Thresholds to simulate
    threshold_results = simulate_thresholds(df, thresholds)

    reports = []

    # Generate the initial report with the new figures
    current_dir = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(current_dir, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template('initial_report_template.html')

    # Render initial report
    initial_html_content = template.render(
        summary_stats=summary_stats,
        umi_read_counts_image=umi_read_counts_image,
        umi_feature_counts_image=umi_feature_counts_image,
        cell_read_counts_image=cell_read_counts_image,
        cell_feature_counts_image=cell_feature_counts_image,
        top_two_ratio_image=top_two_ratio_image
    )

    reports.append(initial_html_content)

    for threshold, df_combined in threshold_results.items():
        print(f"Generating plots for UMI proportion filter threshold {threshold}")

        # Categorize each read-mate after thresholding
        threshold_category_counts = df_combined['threshold_category'].value_counts().reset_index()
        threshold_category_counts.columns = ['category', 'count']

        # Generate bar plot for thresholding categories
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(data=threshold_category_counts, x='category', y='count', ax=ax)
        ax.set_title(f'Thresholding Step: Category Counts for Threshold {threshold}')
        ax.set_xlabel('Category')
        ax.set_ylabel('Number of Read-Mates')
        threshold_category_counts_image = generate_base64_image(fig)
        plt.close(fig)

        # Filter for transitions, exclude eliminations
        df_threshold_filtered = df_combined[
            (df_combined['threshold_category'] == 'reduced features') &
            (df_combined['filtered_features'] != '')
        ]

        # Group on 'original_features' and 'filtered_features', get counts
        threshold_feature_transitions = df_threshold_filtered.groupby(
            ['original_features', 'filtered_features']
        ).size().reset_index(name='count')
        threshold_feature_transitions = threshold_feature_transitions.sort_values('count', ascending=False)

        # Get top N transitions
        N = 20
        top_threshold_transitions = threshold_feature_transitions.head(N)

        # Get top N original features that resulted in zero passing hits after thresholding
        zero_passing_threshold = df_combined[df_combined['threshold_category'] == 'zero passing threshold']
        zero_passing_threshold_features = zero_passing_threshold['original_features'].value_counts().reset_index()
        zero_passing_threshold_features.columns = ['original_features', 'count']
        top_zero_passing_threshold_features = zero_passing_threshold_features.head(N)

        # For intersection step, counts are at UMI level
        # Need to get unique UMIs and their intersection categories
        intersection_category_counts = df_combined[['cb', 'umi', 'intersection_category']].drop_duplicates()
        intersection_category_counts = intersection_category_counts['intersection_category'].value_counts().reset_index()
        intersection_category_counts.columns = ['category', 'count']

        # Generate bar plot for intersection categories
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(data=intersection_category_counts, x='category', y='count', ax=ax)
        ax.set_title(f'Intersection Step: Category Counts for Threshold {threshold}')
        ax.set_xlabel('Category')
        ax.set_ylabel('Number of UMIs')
        intersection_category_counts_image = generate_base64_image(fig)
        plt.close(fig)

        # Filter for transitions, exclude eliminations
        df_intersection_filtered = df_combined[
            (df_combined['intersection_category'] == 'reduced features') &
            (df_combined['post_intersection_features'] != '')
        ].drop_duplicates(subset=['cb', 'umi'])

        # Group on 'original_features' and 'post_intersection_features', get counts
        intersection_feature_transitions = df_intersection_filtered.groupby(
            ['original_features', 'post_intersection_features']
        ).size().reset_index(name='count')
        intersection_feature_transitions = intersection_feature_transitions.sort_values('count', ascending=False)

        # Get top N transitions
        top_intersection_transitions = intersection_feature_transitions.head(N)

        # Get top N original features that resulted in zero passing hits after intersection
        zero_passing_intersection = df_combined[df_combined['intersection_category'] == 'zero passing intersection']
        zero_passing_intersection = zero_passing_intersection[['cb', 'umi', 'original_features']].drop_duplicates()
        zero_passing_features = zero_passing_intersection['original_features'].value_counts().reset_index()
        zero_passing_features.columns = ['original_features', 'count']
        top_zero_passing_features = zero_passing_features.head(N)

        # Load Jinja2 template
        template = env.get_template('threshold_report_template.html')

        # Render report
        html_content = template.render(
            threshold=threshold,
            threshold_category_counts_image=threshold_category_counts_image,
            top_threshold_transitions=top_threshold_transitions.to_dict(orient='records'),
            top_zero_passing_threshold_features=top_zero_passing_threshold_features.to_dict(orient='records'),
            intersection_category_counts_image=intersection_category_counts_image,
            top_intersection_transitions=top_intersection_transitions.to_dict(orient='records'),
            top_zero_passing_features=top_zero_passing_features.to_dict(orient='records')
        )

        reports.append(html_content)

    return reports

def simulate_thresholds(df, thresholds):
    results = {}
    for threshold in thresholds:
        print(f"Simulating UMI count proportion threshold: {threshold}")
        df_copy = df.copy()

        if threshold == 0:
            # Disable thresholding
            df_copy['filtered_features'] = df_copy['features']
        else:
            df_copy = per_umi_thresholding(df_copy, threshold)

        # Collect features after thresholding at the read-mate level
        df_thresholded = df_copy[['cb', 'umi', 'features', 'filtered_features']].copy()

        # Merge original and post-threshold features
        df_combined_threshold = pd.merge(
            df[['cb', 'umi', 'features']].rename(columns={'features': 'original_features'}),
            df_thresholded[['cb', 'umi', 'features', 'filtered_features']],
            left_on=['cb', 'umi', 'original_features'],
            right_on=['cb', 'umi', 'features'],
            how='left'
        )

        # Drop redundant 'features' column
        df_combined_threshold.drop(columns=['features'], inplace=True)

        # Fill NaNs with empty strings
        df_combined_threshold['filtered_features'] = df_combined_threshold['filtered_features'].fillna('')

        # Ensure features are sorted
        df_combined_threshold['original_features'] = df_combined_threshold['original_features'].apply(
            lambda x: ','.join(sorted(x.split(','))) if x else ''
        )
        df_combined_threshold['filtered_features'] = df_combined_threshold['filtered_features'].apply(
            lambda x: ','.join(sorted(x.split(','))) if x else ''
        )

        # Categorize each read-mate after thresholding
        def categorize_threshold(row):
            if row['filtered_features'] == '':
                return 'zero passing threshold'
            elif row['original_features'] == row['filtered_features']:
                return 'unchanged'
            else:
                return 'reduced features'

        df_combined_threshold['threshold_category'] = df_combined_threshold.apply(categorize_threshold, axis=1)

        df_intersected = umi_intersection(df_copy)
        df_post = df_intersected[['cb', 'umi', 'filtered_features']].copy()
        df_post.rename(columns={'filtered_features': 'post_intersection_features'}, inplace=True)

        df_combined = pd.merge(
            df_combined_threshold,
            df_post,
            on=['cb', 'umi'],
            how='left'
        )

        # Fill NaNs with empty strings
        df_combined['post_intersection_features'] = df_combined['post_intersection_features'].fillna('')

        # Ensure features are sorted
        df_combined['post_intersection_features'] = df_combined['post_intersection_features'].apply(
            lambda x: ','.join(sorted(x)) if isinstance(x, list) else (','.join(sorted(x.split(','))) if x else '')
        )

        # Categorize each read-mate after intersection
        def categorize_intersection(row):
            if row['post_intersection_features'] == '':
                return 'zero passing intersection'
            elif row['original_features'] == row['post_intersection_features']:
                return 'unchanged'
            else:
                return 'reduced features'

        df_combined['intersection_category'] = df_combined.apply(categorize_intersection, axis=1)

        results[threshold] = df_combined

    return results

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
