#!/usr/bin/env python
import sys
import os
import shutil
import ctypes as ct
import csv

import argparse
import subprocess
import requests
import stat
import json
import pathlib
import zipfile
import pysam
import pandas as pd
import numpy as np
from collections import Counter
from functools import reduce
from importlib.metadata import version as get_version 
import gzip

from sys import platform

from nimble.types import Config
from nimble.parse import parse_fasta, parse_filter_config, parse_csv
from nimble.usage import print_usage_and_exit
from nimble.utils import get_exec_name_from_platform, append_path_string, per_umi_thresholding, umi_intersection
from nimble.report_generation import generate_plots

ALIGN_TRIES = 10
ALIGN_TRIES_THRESHOLD = 0

def validate_gzip(file_path):
    try:
        with gzip.open(file_path, 'rb') as f:
            while f.read(1024):  # Read in chunks to validate
                pass
        print(f"GZIP validation successful for {file_path}")
    except Exception as e:
        print(f"GZIP validation failed for {file_path}: {e}")

# Generate and write human-editable config json files to disk. Input data is a CSV, FASTA, or both
def generate(file, opt_file, output_path):

    (data, config, is_csv_req) = process_file(file, opt_file)
    (data_opt, config_opt, is_csv_opt) = process_file(opt_file, file)

    final_config = config
    if data_opt != None and is_csv_opt:
        final_config = config_opt

    final_data = None
    if data_opt != None:
        if is_csv_req:
            final_data = collate_data(data_opt, data)
        elif is_csv_opt:
            final_data = collate_data(data, data_opt)
    else:
        final_data = data

    # Write reference and default config to disk
    with open(output_path, "w") as f:
        json.dump([ final_config.__dict__, final_data.__dict__], f, indent=2)

# Parse a lone FASTA/lone CSV/CSV-FASTA tuple. If there's a lone FASTA, generate a simple library.
# If there's a lone CSV, assume it has sequence information or a genbank link.
# If there is not a lone CSV, assume it contains the metadata and that the sequence information is contained in the FASTA.
def process_file(file, paired_file):
    data = None
    config = None
    is_csv = False

    if file:
        if pathlib.Path(file).suffix == ".fasta":
            (data, config) = parse_fasta(file)
        elif pathlib.Path(file).suffix == ".csv" and paired_file:
            (data, config) = parse_csv(file, False)
            is_csv = True
        elif pathlib.Path(file).suffix == ".csv" and not paired_file:
            (data, config) = parse_csv(file, True)
            is_csv = True

    return (data, config, is_csv)


def collate_data(data, metadata):
    name_idx = data.headers.index("sequence_name")
    sequence_idx = data.headers.index("sequence")
    nt_length_idx = data.headers.index("nt_length")

    meta_name_idx = metadata.headers.index("sequence_name")
    meta_sequence_idx = metadata.headers.index("sequence")
    meta_nt_length_idx = metadata.headers.index("nt_length")

    metadata.columns[meta_sequence_idx] = ["" for _ in range(0, len(data.columns[sequence_idx]))]
    metadata.columns[meta_nt_length_idx] = ["" for _ in range(0, len(data.columns[nt_length_idx]))]

    for (from_idx, name) in enumerate(data.columns[name_idx]):
        if name not in metadata.columns[meta_name_idx]:
            print("Error -- record " + name + " is not found in both input files.")
            sys.exit()

        update_idx = metadata.columns[meta_name_idx].index(name)

        metadata.columns[meta_sequence_idx][update_idx] = data.columns[sequence_idx][from_idx]
        metadata.columns[meta_nt_length_idx][update_idx] = data.columns[nt_length_idx][from_idx]

    return metadata


# Given the name of a release, download the platform-specific executable from that release.
# Given no name, default to the most recent release.
def download(release):
    exec_name = get_exec_name_from_platform()
    print("Downloading " + exec_name)

    url = ""
    if len(release) == 1:
        url = (
            "https://github.com/BimberLab/nimble-aligner/releases/download/"
            + release[0]
            + "/"
            + exec_name
        )
    else:
        url = (
            "https://github.com/BimberLab/nimble-aligner/releases/latest/download/"
            + exec_name
        )

    # Download aligner
    r = requests.get(url)

    # The executable will be placed in the python package directory
    aligner_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")
    print("Aligner download path: " + aligner_path)

    # If we successfully downloaded it, write the file to disk and give it execution permissions if necessary
    if r.status_code == 200:
        with open(aligner_path, "wb") as f:
            f.write(r.content)
        if sys.platform == "linux" or sys.platform == "darwin": # If we're on a Unix, ensure permission to write
            st = os.stat(aligner_path)
            os.chmod(aligner_path, st.st_mode | stat.S_IEXEC | stat.S_IXOTH)
    else:
        print("Error -- could not download aligner, status code " + str(r.status_code))
        sys.exit()


# Check if the aligner exists -- if it does, call it with the given parameters.
def align(reference, output, input, num_cores, strand_filter, trim, tmpdir):
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")

    if not os.path.exists(path):
        print("No aligner found. Attempting to download the latest release.\n")
        download([])

        global ALIGN_TRIES
        ALIGN_TRIES = ALIGN_TRIES + 1
        if ALIGN_TRIES >= ALIGN_THRESHOLD:
            print("Error -- could not find or download aligner.")
            sys.exit()

        return align(reference, output, input, num_cores, strand_filter, trim, tmpdir)

    print("Aligning input data to the reference libraries")
    sys.stdout.flush()

    input_ext = os.path.splitext(input[0])[-1].lower()
    should_delete_input = False
    if input_ext == ".bam":
        input[0] = sort_input_bam(input[0], num_cores, tmpdir)
        should_delete_input = True

    processed_param_list = []
    for input_file in input:
        processed_param_list.extend(["--input", input_file])
    processed_param_list.extend(["-c", str(num_cores), "--strand_filter", strand_filter])

    library_list = reference.split(",")
    for library in library_list:
        out_file_append = ""

        if len(library_list) > 1:
            out_file_append = "." + os.path.splitext(os.path.basename(library))[0]

        processed_param_list.extend(["-r", library, "-o", append_path_string(output, out_file_append)])

    if trim != "":
        processed_param_list.extend(["-t", trim])

    print(processed_param_list)

    proc = subprocess.Popen([path] + processed_param_list)
    proc.wait()

    return_code = proc.returncode

    if should_delete_input and return_code == 0:
        print("Deleting intermediate sorted .bam file")

        try:
            os.remove(input[0])
            os.remove(input[0] + ".done")
        except Exception as e:
            print(f"Error when attempting to delete sorted .bam file: {e}")
    else:
        print("Retaining all input files.")

    return return_code

def report(input, output, summarize_columns_list=None, threshold=0.05, disable_thresholding=False):
    df = None
    
    # If the file has data, try to read it. Write an empty output and return if there is no data.
    if os.path.getsize(input) > 0:
        try:
            df = pd.read_csv(input, sep='\t', quoting=csv.QUOTE_NONE, engine='python')
            # Use the r1 version of the cb and umi flags
            df.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi', 'nimble_features': 'features'}, inplace=True)
        except pd.errors.EmptyDataError:
            write_empty_df(output)
            return
    else:
        write_empty_df(output)
        return

    # If the file is not empty but the DataFrame is, write an empty output
    if df.empty:
        write_empty_df(output)
        return

    # Keep only necessary columns
    df_init = df.copy()
    df = df[['features', 'umi', 'cb', 'nimble_score']]

    # Drop rows where 'features', 'umi', 'cb', or 'nimble_score' are null or empty
    df = df.dropna(subset=['features', 'umi', 'cb', 'nimble_score'])
    df = df[(df['features'] != '') & (df['umi'] != '') & (df['cb'] != '')]

    # Ensure features are sorted
    df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))

    # Merge duplicate rows after sorting ambiguous features, summing 'nimble_score'
    df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()

    # If the file is not empty but the DataFrame is after filtering, write an empty output
    if df.empty:
        write_empty_df(output)
        return

    # Apply thresholding algorithm unless disabled
    if not disable_thresholding:
        df = per_umi_thresholding(df, threshold)
    else:
        df['filtered_features'] = df['features']

    # Apply UMI intersection
    df_grouped = umi_intersection(df)

    # Identify and drop rows with empty intersections
    empty_intersection_rows = df_grouped['filtered_features'].apply(lambda x: len(x) == 0)
    empty_intersection_count = empty_intersection_rows.sum()
    print(f"Dropped {empty_intersection_count} UMIs due to empty intersections")
    df_grouped = df_grouped[~empty_intersection_rows]

    # Convert intersected features back to strings
    df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))

    # Rename columns
    df_grouped.columns = ['cell_barcode', 'umi', 'feature']

    # For the final counts, we count the number of UMIs per cell barcode and feature
    df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
    df_counts = df_counts[['feature', 'count', 'cell_barcode']]

    # Write to output file
    df_counts.to_csv(output, sep='\t', index=False, header=False)

    if summarize_columns_list:
        summary_output = "summarize." + output
        summarize_fields(df_init, summarize_columns_list, summary_output)

def write_empty_df(output):
    print('No data to parse from input file, writing empty output.')
    empty_df = pd.DataFrame(columns=['feature', 'count', 'cell_barcode'])
    empty_df.to_csv(output, sep='\t', index=False, header=False)

def summarize_fields(df, columns, output_file):
    summary_df = df.groupby('umi')[columns].agg(lambda x: x.value_counts().to_dict())
    summary_df = summary_df.applymap(lambda x: '; '.join([f'{k}({v})' for k, v in x.items()]))
    summary_df.reset_index().to_csv(output_file, sep='\t', index=False)

def sort_input_bam(bam, cores, tmp_dir):
    print("Sorting input bam")

    if not tmp_dir:
        tmp_dir = os.environ.get("TMPDIR")

    created_tmp_dir = False

    if tmp_dir:
        samtools_tmp = os.path.join(tmp_dir, "samtools_tmp")
    else:
        samtools_tmp = None

    sorted_bam = os.path.join(tmp_dir, "sorted-" + os.path.basename(bam))
    sorted_bam_done = sorted_bam + ".done"

    print("Sorting " + bam + " Outputting to " + sorted_bam)
    sys.stdout.flush()

    if os.path.exists(sorted_bam_done):
        print("Sorted bam file already exists, skipping the sorting step.")
        sys.stdout.flush()
        return sorted_bam

    if tmp_dir and not os.path.exists(samtools_tmp):
        try:
            os.makedirs(samtools_tmp)
            created_tmp_dir = True
            print(f"Created temporary directory {samtools_tmp}")
        except OSError as e:
            print(f"Could not create temporary directory {samtools_tmp}: {e}")
            sys.stdout.flush()

    if tmp_dir:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', str(cores), '-T', samtools_tmp, bam)
    else:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', str(cores), bam)

    sort_log = pysam.sort.get_messages()

    if (sort_log):
        print("samtools messages: " + sort_log)
        
    if created_tmp_dir:
        try:
            shutil.rmtree(samtools_tmp)
            print(f"Deleted samtools temporary directory {samtools_tmp}")
        except Exception as e:
            print(f"Could not delete temporary directory {samtools_tmp}: {e}")
        sys.stdout.flush()

    # Note: this allows the code to know if sorting completed successfully, as opposed to dying in the middle
    with open(sorted_bam_done, 'w') as fp:
        pass

    if not os.path.exists(sorted_bam_done):
        raise RuntimeError("unable to create BAM done_file")

    return sorted_bam


if __name__ == "__main__":
    csv.field_size_limit(int(ct.c_ulong(-1).value // 2))

    try:
        nimble_version = get_version('nimble')
    except Exception:
        nimble_version = 'unknown'

    parser = argparse.ArgumentParser(description='nimble align')
    parser.add_argument('-v', '--version', action='version', version=f'nimble {nimble_version}')

    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    download_parser = subparsers.add_parser('download')
    download_parser.add_argument('--release', help='The release to download.', type=str, default=[])

    generate_parser = subparsers.add_parser('generate')
    generate_parser.add_argument('--file', help='The file to process.', type=str, required=True)
    generate_parser.add_argument('--opt-file', help='The optional file to process.', type=str, default=None)
    generate_parser.add_argument('--output_path', help='The path to the output file.', type=str, required=True)

    align_parser = subparsers.add_parser('align')
    align_parser.add_argument('--reference', help='The reference genome to align to.', type=str, required=True)
    align_parser.add_argument('--output', help='The path to the output file.', type=str, required=True)
    align_parser.add_argument('--input', help='The input reads.', type=str, required=True, nargs='+')
    align_parser.add_argument('-c', '--num_cores', help='The number of cores to use for alignment.', type=int, default=1)
    align_parser.add_argument('--strand_filter', help='Filter reads based on strand information.', type=str, default="unstranded")
    align_parser.add_argument('--trim', help='Configuration for trimming read-data, in the format <TARGET_LENGTH>:<STRICTNESS>, comma-separated, one entry for each passed library', type=str, default="")
    align_parser.add_argument('--tmpdir', help='Path to a temporary directory for sorting .bam files', type=str, default=None)

    report_parser = subparsers.add_parser('report')
    report_parser.add_argument('-i', '--input', help='The input file.', type=str, required=True)
    report_parser.add_argument('-o', '--output', help='The path to the output file.', type=str, required=True)
    report_parser.add_argument('-s', '--summarize', help='CSV list of columns to summarize.', type=str, default=None)
    report_parser.add_argument(
        '-t', '--threshold', 
        help='Proportional count threshold for filtering features (default: 0.05).', 
        type=float, 
        default=0.05
    )
    report_parser.add_argument(
        '--disable_thresholding',
        help='Disable the per-UMI proportional count thresholding algorithm.',
        action='store_true',
        default=False
    )

    plot_parser = subparsers.add_parser('plot')
    plot_parser.add_argument('--input_file', help='The nimble counts output file to process.', type=str, required=True)
    plot_parser.add_argument('--output_file', help='Filepath for the HTML report.', type=str, required=True)

    args = parser.parse_args()

    if args.subcommand == 'download':
        download(args.release)
    elif args.subcommand == 'generate':
        generate(args.file, args.opt_file, args.output_path)
    elif args.subcommand == 'align':
        sys.exit(align(args.reference, args.output, args.input, args.num_cores, args.strand_filter, args.trim, args.tmpdir))
    elif args.subcommand == 'report':
        summarize_columns_list = args.summarize.split(',') if args.summarize else None
        report(args.input, args.output, summarize_columns_list, args.threshold, args.disable_thresholding)

    elif args.subcommand == 'plot':
        if os.path.getsize(args.input_file) > 0:
            try:
                print(f"Loading alignment data from {args.input_file}")
                df = pd.read_csv(args.input_file, sep='\t', quoting=csv.QUOTE_NONE, engine='python')
                generate_plots(df, args.output_file)
            except pd.errors.EmptyDataError:
                print("Input file is empty.")
        else:
            print("Input file is empty.")
