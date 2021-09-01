#!/usr/bin/env python
import sys
import os
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

align_tries = 10
align_tries_threshold = 0


# Generate and write human-editable config json files to disk
def generate(file, output_path):
    #TODO process FASTA + CSV
    data = None
    config = None

    if pathlib.Path(file).suffix == ".fasta":
        (data, config) = parse_fasta(file)
    elif pathlib.Path(file).suffix == ".csv":
        (data, config) = parse_csv(file)

    # Write reference and default config to disk
    with open(output_path, "w") as f:
        json.dump([ config.__dict__, data.__dict__], f, indent=2)


# Given the name of a release, download the platform-specific executable from that release.
# Given no name, default to the most recent release.
def download(release):
    exec_name = get_exec_name_from_platform()

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
            os.chmod(aligner_path, st.st_mode | stat.S_IEXEC)
    else:
        print("Error -- could not download aligner, status code " + str(r.status_code))
        sys.exit()


# Check if the aligner exists -- if it does, call it with the given parameters.
def align(param_list):
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")
    input_ext = os.path.splitext(param_list[2])[-1].lower()

    cores = "1"
    if "-c" in param_list:
        cores = param_list[param_list.index("-c") + 1]
    elif "--cores" in param_list:
        cores = param_list[param_list.index("--cores") + 1]

    if os.path.exists(path):
        if input_ext == ".bam":
            sort_input_bam(param_list[2], cores)
        print("Aligning input .bam to the reference library")
        subprocess.call([path] + param_list)
    else:
        print("No aligner found. Attempting to download the latest release.\n")
        download_aligner([],)

        align_tries += 1
        if align_tries >= align_threshold:
            print("Error -- could not find or download aligner.")
            sys.exit()

        align(param_list)


def sort_input_bam(bam, cores):
    print("Sorting input .bam")
    tmp_dir = os.environ.get("TMPDIR")

    if tmp_dir:
        pysam.sort('-t', 'UR', '-n', '-o', bam, '-@', cores, bam, '-T', tmp_dir)
    else:
        pysam.sort('-t', 'UR', '-n', '-o', bam, '-@', cores, bam)


if __name__ == "__main__":
    if len(sys.argv) == 1:  # Ensure we can index sys.argv[1]
        print_usage_and_exit()
    elif sys.argv[1] == "download" and len(sys.argv) <= 3:
        download(sys.argv[2:])
    elif sys.argv[1] == "generate" and len(sys.argv) == 4:
        generate(sys.argv[2], sys.argv[3])
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
