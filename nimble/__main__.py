#!/usr/bin/env python
import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
import subprocess
import requests
import stat
import json
import pathlib
import zipfile
import pysam

from sys import platform

from nimble.types import Config
from nimble.parse import parse_fasta, parse_filter_config, parse_csv
from nimble.usage import print_usage_and_exit
from nimble.utils import get_exec_name_from_platform
from nimble.reporting import report

ALIGN_TRIES = 10
ALIGN_TRIES_THRESHOLD = 0


# Generate and write human-editable config json files to disk
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

    aligner_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")

    # If we successfully downloaded it, write the file to disk and give it execution permissions if necessary
    if r.status_code == 200:
        with open(aligner_path, "wb") as f:
            f.write(r.content)
        if sys.platform == "linux" or sys.platform == "darwin":
            st = os.stat(aligner_path)
            os.chmod(aligner_path, st.st_mode | stat.S_IEXEC | stat.S_IXOTH)
    else:
        print("Error -- could not download aligner, status code " + str(r.status_code))
        sys.exit()


# Check if the aligner exists -- if it does, call it with the given parameters.
def align(param_list):
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")

    param_list_noflags = []
    it = iter(param_list)
    for param in it:
        if param[0] == "-":
            next(it, None)
        else:
            param_list_noflags.append(param)

    bam_param_idx = param_list.index(param_list_noflags[2])
    input_ext = os.path.splitext(param_list_noflags[2])[-1].lower()

    cores = "1"
    if "-c" in param_list:
        cores = param_list[param_list.index("-c") + 1]
    elif "--cores" in param_list:
        cores = param_list[param_list.index("--cores") + 1]

    if os.path.exists(path):
        if input_ext == ".bam":
            split = param_list_noflags[2].rsplit("/", 1)
            sort_input_bam(split, cores)
            param_list[bam_param_idx] = split[0] + "/sorted-" + split[1]

        print("Aligning input .bam to the reference library")
        sys.stdout.flush()
        subprocess.call([path] + param_list)

        print("Deleting intermediate sorted .bam file")
        os.remove(param_list[bam_param_idx])
    else:
        print("No aligner found. Attempting to download the latest release.\n")
        download([])

        global ALIGN_TRIES
        ALIGN_TRIES = ALIGN_TRIES + 1
        if ALIGN_TRIES >= ALIGN_THRESHOLD:
            print("Error -- could not find or download aligner.")
            sys.exit()

        align(param_list)


def sort_input_bam(file_tuple, cores):
    print("Sorting input .bam")
    sys.stdout.flush()
    tmp_dir = os.environ.get("TMPDIR")

    bam = file_tuple[0] + "/" + file_tuple[1]
    sorted_bam = file_tuple[0] + "/sorted-" + file_tuple[1]

    if tmp_dir:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', cores, bam, '-T', tmp_dir)
    else:
        pysam.sort('-t', 'UR', '-n', '-o', sorted_bam, '-@', cores, bam)


if __name__ == "__main__":
    if len(sys.argv) == 1:  # Ensure we can index sys.argv[1]
        print_usage_and_exit()
    elif sys.argv[1] == "download" and len(sys.argv) <= 3:
        download(sys.argv[2:])
    elif sys.argv[1] == "generate" and len(sys.argv) >= 4 and len(sys.argv) <= 5:
        if len(sys.argv) == 4:
            generate(sys.argv[2], None, sys.argv[3])
        else:
            generate(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "align":
        align(sys.argv[2:])
    elif sys.argv[1] == "report" and len(sys.argv) == 5:
        (methods, values) = parse_filter_config(sys.argv[2])
        report(methods, values, sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "filter" and len(sys.argv) >= 5 and len(sys.argv) <= 6:
        if len(sys.argv) == 5:
            report([sys.argv[2]], [None], sys.argv[3], sys.argv[4])
        else:
            report([sys.argv[2]], [int(sys.argv[3])], sys.argv[4], sys.argv[5])
    else:
        print_usage_and_exit()
