import sys
import os
import distro
import pandas as pd
import numpy as np

from pathlib import Path

def append_path_string(input_path, pathAppendString):
    # Extract the filename from the input path
    filename = os.path.basename(input_path)

    # Split the filename into the root and extension
    root = filename
    ext = ""
    while True:
        root, ext2 = os.path.splitext(root)
        if ext2 == "":
            break
        ext = ext2 + ext

    # Concatenate the root and the pathAppendString
    new_filename = root + pathAppendString + ext

    # Join the new filename with the path to form the full path
    new_path = os.path.join(os.path.dirname(input_path), new_filename)
    return new_path


# Generate reference genome name by getting the filename and prettifying it
def get_library_name_from_filename(seq_path):
    return Path(seq_path).stem.replace("_", " ")


# Get the name of a Github release given the target platform
def get_exec_name_from_platform():
    exec_name = ""

    if sys.platform == "win32":
        exec_name = "Windows.exe"
    elif sys.platform == "linux":
        exec_name = "Linux.out"
    elif sys.platform == "darwin":
        exec_name = "MacOS.out"
    else:
        print(
            "Error -- platform "
            + sys.platform
            + " does not have a compatible release on GitHub. Please open an issue at https://github.com/BimberLab/nimble."
        )
        sys.exit()

    return exec_name


# Serialize results of an alignment/filter to output TSV, while collapsing equal ambiguity classes
def write_data_to_tsv(output_path, out_data, metadata):
    def _modify_metadata_score(metadata, new_score):
        metadata[0] = str(new_score)

        if len(current_metadata) == 1:
            metadata[0] += "\n"

    def _add_line_to_cache(cache, row, current_metadata, current_score):
        if len(row) > 0:
            key = ",".join([reference for reference in row])

            if key in cache:
                cache_score = int(cache[key].split("\t")[0])
                cache_score += current_score
                _modify_metadata_score(current_metadata, cache_score)

            line_metadata = "\t".join([elem for elem in current_metadata])
            cache[key] = line_metadata

    cache = {}
    metadata = iter(metadata)
    scores = iter(out_data[0])

    out_data = out_data.drop(0, axis=1)
    out_data = iter(
        out_data.apply(lambda row: row.values[~pd.isna(row.values)], axis=1)
    )

    str_rep = ""
    current_row = next(out_data)
    current_score = next(scores)
    current_metadata = None

    with open(output_path, "w") as f:
        f.write(next(metadata))

        if current_metadata == None:
            current_metadata = next(metadata)

        for row in out_data:
            if not np.array_equal(current_row, row):
                _modify_metadata_score(current_metadata, current_score)
                _add_line_to_cache(cache, current_row, current_metadata, current_score)

                current_row = row
                current_score = next(scores)
                current_metadata = next(metadata)
                continue

            current_score += next(scores)
            current_metadata = next(metadata)

        _add_line_to_cache(cache, current_row, current_metadata, current_score)

        for key, value in cache.items():
            str_rep += key + "\t" + value

        f.write(str_rep)

def trim_low_complexity_regions(seq):
    return seq

def per_umi_thresholding(df, threshold):
    def filter_umi_features(umi_group):
        counts = []
        total_score = 0

        # Initial proportional scores
        for _, row in umi_group.iterrows():
            features = row['features'].split(',')
            score_per_feature = row['nimble_score'] / len(features)
            total_score += row['nimble_score']
            for feature in features:
                counts.append({'feature': feature, 'nimble_score': score_per_feature})

        # Aggregate scores per feature
        feature_scores = pd.DataFrame(counts).groupby('feature')['nimble_score'].sum()

        # Initialize filtered_features_set
        filtered_features_set = None

        # Iterative filtering
        while True:
            # Check to handle empty feature_scores
            if feature_scores.empty:
                # No features left after filtering
                filtered_features_set = set()
                break

            feature_ratios = feature_scores / total_score
            to_drop = feature_ratios[feature_ratios < threshold].index

            if len(to_drop) == 0:
                # No more features to drop
                filtered_features_set = set(feature_scores.index)
                break

            # Reassign scores, excluding filtered features
            filtered_counts = []
            total_score = 0

            for _, row in umi_group.iterrows():
                features = [f for f in row['features'].split(',') if f not in to_drop]
                if not features:
                    continue
                score_per_feature = row['nimble_score'] / len(features)
                total_score += row['nimble_score']
                for feature in features:
                    filtered_counts.append({'feature': feature, 'nimble_score': score_per_feature})

            # Check to handle empty filtered_counts
            if not filtered_counts:
                # No features left after filtering
                filtered_features_set = set()
                break

            feature_scores = pd.DataFrame(filtered_counts).groupby('feature')['nimble_score'].sum()

        # Prepare the result DataFrame for this UMI group
        result = umi_group.copy()
        result['filtered_features'] = ''

        # Ensure filtered_features_set is not None
        if filtered_features_set is None:
            filtered_features_set = set()

        # Assign filtered features per read-mate
        def assign_filtered_features(row):
            original_features = set(row['features'].split(','))
            features_to_keep = original_features & filtered_features_set
            return ','.join(sorted(features_to_keep)) if features_to_keep else ''

        result['filtered_features'] = result.apply(assign_filtered_features, axis=1)

        return result[['cb', 'umi', 'features', 'filtered_features']]

    # Apply the function to each UMI group
    filtered_umis = df.groupby(['cb', 'umi'], group_keys=False).apply(filter_umi_features).reset_index(drop=True)

    # Merge the filtered results back to the main DataFrame on 'cb', 'umi', 'features'
    df = pd.merge(
        df,
        filtered_umis[['cb', 'umi', 'features', 'filtered_features']],
        on=['cb', 'umi', 'features'],
        how='inner'
    )

    # Filter out rows with empty 'filtered_features'
    df = df[df['filtered_features'] != '']

    return df

def umi_intersection(df):
    # Convert filtered_features to lists
    df['filtered_features'] = df['filtered_features'].str.split(',')

    # Group by cb and umi, collect lists of lists
    df_grouped = df.groupby(['cb', 'umi'], group_keys=False)['filtered_features'].apply(list).reset_index()

    # Apply intersection to the lists of lists
    df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(intersect_lists)

    return df_grouped

def intersect_lists(list_of_lists):
    if not list_of_lists:
        return []
    return sorted(set.intersection(*map(set, list_of_lists)))
