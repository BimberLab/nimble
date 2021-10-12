import sys
import distro
import pandas as pd
import numpy as np

from pathlib import Path


# Derived from experimental data
LOW_COMPLEXITY_REGION_LEN = 15


# Generate reference genome name by getting the filename and prettifying it
def get_library_name_from_filename(seq_path):
    return Path(seq_path).stem.replace("_", " ")


# Get the name of a Github release given the target platform
def get_exec_name_from_platform():
    exec_name = ""

    if sys.platform == "win32":
        exec_name = "windows.exe"
    elif sys.platform == "linux":
        exec_name = "Linux.out"
    elif sys.platform == "darwin":
        exec_name = "Mac.app"
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
    current_base = None
    current_region = None
    regions = []
    new_seq = ""

    # Partition the sequence into contiguous regions
    for base in seq:
        if current_base == None:
            current_base = base
            current_region = current_base
            continue

        if base != current_base:
            regions.append(current_region)
            current_base = base
            current_region = current_base
        else:
            current_region += base

    regions.append(current_region)

    print(regions)

    # Concat all of the regions, skipping contiguous regions with length >= LOW_COMPLEXITY_REGION_LEN
    for region in regions:
        if len(region) < LOW_COMPLEXITY_REGION_LEN:
            new_seq += region
        else:
            print(region)

    return new_seq
