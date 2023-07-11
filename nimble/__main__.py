#!/usr/bin/env python
import sys
import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'

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

from sys import platform

from nimble.types import Config
from nimble.parse import parse_fasta, parse_filter_config, parse_csv
from nimble.usage import print_usage_and_exit
from nimble.utils import get_exec_name_from_platform, low_complexity_filter_amount, append_path_string

ALIGN_TRIES = 10
ALIGN_TRIES_THRESHOLD = 0


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

    print("Filtered " + str(low_complexity_filter_amount) + " base pairs from reference library.")

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
def align(reference, output, input, _alignment_path, _log_path, num_cores, strand_filter):
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")
    input_ext = os.path.splitext(input)[-1].lower()

    if os.path.exists(path):
        if input_ext == ".bam":
            split = input.rsplit("/", 1)
            sort_input_bam(split, num_cores)
            input = split[0] + "/sorted-" + split[1]


        print("Aligning input .bam to the reference libraries")
        sys.stdout.flush()

        library_list = reference.split(",")

        processed_param_list = ["--input", input, "-c", str(num_cores), "--strand_filter", strand_filter]

        for library in library_list:
            out_file_append = ""

            if len(library_list) > 1:
                out_file_append = "." + os.path.splitext(os.path.basename(library))[0]

            processed_param_list.extend(["-r", library, "-o", append_path_string(output, out_file_append)])
        print(processed_param_list)

        proc = subprocess.Popen([path] + processed_param_list)
        proc.wait()

        return_code = proc.returncode

        if input_ext == ".bam" and return_code == 0:
            print("Deleting intermediate sorted .bam file")
            os.remove(input)

        return return_code
    else:
        print("No aligner found. Attempting to download the latest release.\n")
        download([])

        global ALIGN_TRIES
        ALIGN_TRIES = ALIGN_TRIES + 1
        if ALIGN_TRIES >= ALIGN_THRESHOLD:
            print("Error -- could not find or download aligner.")
            sys.exit()

        return align(reference, output, input, alignment_path, log_path, num_cores, strand_filter)

def intersect_lists(lists):
    # Returns the intersection of all lists in a list
    return list(reduce(set.intersection, map(set, lists)))

def report(input, output):
    # Read input file
    df = pd.read_csv(input, sep='\t', compression='gzip')

    # Keep only necessary columns
    df = df[['features', 'umi', 'cb']]

    # Drop rows where 'features', 'umi', or 'cb' are null or empty
    df = df.dropna(subset=['features', 'umi', 'cb'])
    df = df[(df['features'] != '') & (df['umi'] != '') & (df['cb'] != '')]

    # Split the feature strings into lists
    df['features'] = df['features'].str.split(',')

    # Group by cell barcode (CB) and UMI, aggregate features into a list
    df_grouped = df.groupby(['cb', 'umi'])['features'].apply(list)

    # Calculate the intersection of features within each UMI
    df_grouped = df_grouped.apply(intersect_lists)

    # Convert back to a DataFrame
    df_grouped = df_grouped.reset_index()

    # Identify rows where the intersection resulted in an empty list
    empty_intersection_rows = df_grouped['features'].apply(lambda x: len(x) == 0)

    # Count these rows and print the number
    empty_intersection_count = empty_intersection_rows.sum()
    print(f"Dropped {empty_intersection_count} counts due to empty intersections")

    # Drop these rows from the DataFrame
    df_grouped = df_grouped[~empty_intersection_rows]

    # Join the intersected features back into a string
    df_grouped['features'] = df_grouped['features'].apply(lambda x: ','.join(x))

    # Rename columns
    df_grouped.columns = ['cell_barcode', 'umi', 'feature']

    # Count unique UMIs per cell_barcode-feature pair
    df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index()

    # Rename count column
    df_counts.columns = ['cell_barcode', 'feature', 'count']

    # Reorder the columns
    df_counts = df_counts.reindex(['feature', 'count', 'cell_barcode'], axis=1)

    # Write to output file
    df_counts.to_csv(output, sep='\t', index=False, header=False)

def sort_input_bam(file_tuple, cores):
    print("Sorting input .bam")
    tmp_dir = os.environ.get("TMPDIR")

    bam = ""
    sorted_bam = ""

    if len(file_tuple) > 1:
        bam = file_tuple[0] + "/" + file_tuple[1]
        sorted_bam = file_tuple[0] + "/sorted-" + file_tuple[1]
    else:
        bam = file_tuple[0]
        sorted_bam = "./sorted" + file_tuple[0]

    print("Sorting " + bam + " Outputting to " + sorted_bam)
    sys.stdout.flush()

    if os.path.isfile(sorted_bam):
        print("Sorted bam file already exists, skipping the sorting step.")
        sys.stdout.flush()
        return

    if tmp_dir:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', str(cores), '-T', tmp_dir, bam)
    else:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', str(cores), bam)

    sort_log = pysam.sort.get_messages()

    if (sort_log):
        print("samtools messages: " + sort_log)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
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
    align_parser.add_argument('--input', help='The input reads.', type=str, required=True)
    align_parser.add_argument('--alignment_path', help='The path to the alignment file.', type=str, default=None)
    align_parser.add_argument('--log_path', help='The path to the log file.', type=str, default=None)
    align_parser.add_argument('-c', '--num_cores', help='The number of cores to use for alignment.', type=int, default=1)
    align_parser.add_argument('--strand_filter', help='Filter reads based on strand information.', type=str, default="unstranded")

    report_parser = subparsers.add_parser('report')
    report_parser.add_argument('-i', '--input', help='The input file.', type=str, required=True)
    report_parser.add_argument('-o', '--output', help='The path to the output file.', type=str, required=True)

    args = parser.parse_args()

    if args.subcommand == 'download':
        download(args.release)
    elif args.subcommand == 'generate':
        generate(args.file, args.opt_file, args.output_path)
    elif args.subcommand == 'align':
        sys.exit(align(args.reference, args.output, args.input, args.alignment_path, args.log_path, args.num_cores, args.strand_filter))
    elif args.subcommand == 'report':
        report(args.input, args.output)
