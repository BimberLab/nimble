#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import requests
import stat
import json
import zipfile
import pysam

from sys import platform

from nimble.types import Config
from nimble.parse import parse_fasta, parse_filter_config
from nimble.usage import print_usage_and_exit
from nimble.utils import get_exec_name_from_platform
from nimble.reporting import report

align_tries = 10
align_tries_threshold = 0

# Take human-editable config json files and compile into a single minified file for the aligner
def compile_config(reference_output_path, config_output_path, compiled_json_path):
    with open(reference_output_path, "r") as ref, open(
        config_output_path, "r"
    ) as conf, open(compiled_json_path, "w") as comp:
        reference = json.load(ref)
        config = json.load(conf)
        json.dump([config, reference], comp)


# Generate and write human-editable config json files to disk
def create_editable_config(seq_path, reference_output_path, config_output_path):
    # Determine if we're working with .fasta or a zipped .bam
    file_ext = os.path.splitext(seq_path)[1]

    data = None
    config = None

    (data, config) = parse_fasta(seq_path)

    # Write reference and default config to disk
    with open(reference_output_path, "w") as f:
        json.dump(data.__dict__, f, indent=2)

    with open(config_output_path, "w") as f:
        json.dump(config.__dict__, f, indent=2)


# Given the name of a release, download the platform-specific executable from that release.
# Given no name, default to the most recent release.
def download_aligner(release):
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

    if os.path.exists(path):
        if input_ext == ".bam":
            sort_input_bam(param_list[2])
        print("Aligning input .bam to the reference library")
        subprocess.call([path] + param_list)
    else:
        print("No aligner found. Attempting to download the latest release.\n")
        download_aligner([])

        align_tries += 1
        if align_tries >= align_threshold:
            print("Error -- could not find or download aligner.")
            sys.exit()

        align(param_list)


def sort_input_bam(bam):
    print("Sorting input .bam")
    pysam.sort('-t', 'UR', '-n', '-o', bam, '-@', '8', bam)


if __name__ == "__main__":
    if len(sys.argv) == 1:  # Ensure we can index sys.argv[1]
        print_usage_and_exit()
    elif sys.argv[1] == "download" and len(sys.argv) <= 3:
        download_aligner(sys.argv[2:])
    elif sys.argv[1] == "generate" and len(sys.argv) == 5:
        create_editable_config(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "compile" and len(sys.argv) == 5:
        compile_config(sys.argv[2], sys.argv[3], sys.argv[4])
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
