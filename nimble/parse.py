import json
import csv

import pandas as pd

from Bio import SeqIO
from io import StringIO

from nimble.types import Data, Config, DataType
from nimble.utils import get_library_name_from_filename
from nimble.remote import fetch_sequence, get_ids


# Given a path to a .fasta, reads it and returns a (Data, Config) tuple with filled objects
def parse_fasta(seq_path):
    data = Data()
    config = Config()
    config.data_type = DataType.FASTA

    reference_name = get_library_name_from_filename(seq_path)

    f = SeqIO.parse(seq_path, "fasta")
    for record in f:
        data.columns[0].append(reference_name)

        if record.id is not None:
            data.columns[1].append(record.id)
        else:
            data.columns[1].append("null")

        data.columns[2].append(str(len(record)))
        data.columns[3].append(str(record.seq))

    return (data, config)


# Read data from the backend aligner's output format -- a TSV
def parse_alignment_results(input_path):
    with open(input_path, "r") as f:
        metadata = [next(f)]

        str_rep = ""
        max_line_len = 0

        for line in f:
            csv_line = line.split("\t")[1].strip() + "," + line.split("\t")[0] + "\n"
            str_rep += csv_line + "\n"
            curr_line_len = len(csv_line.split(","))

            if curr_line_len > max_line_len:
                max_line_len = curr_line_len

            metadata.append(line.split("\t")[1:])

    names = [i for i in range(0, max_line_len)]
    return (pd.read_csv(StringIO(str_rep), header=None, names=names), metadata)


# Parse the reference.json format for the list of filters and their configurations
def parse_filter_config(reference_path):
    methods = []
    values = []

    with open(reference_path) as ref:
        data = json.load(ref)

    for method in data[0]["filters"]:
        methods.append(method["name"])
        values.append(method["value"])

    return (methods, values)


def parse_csv(csv_path, has_sequences=True):
    data = Data()
    config = Config()
    config.data_type = DataType.FASTA

    reference_genome = get_library_name_from_filename(csv_path)
    reference_genomes = []
    sequence_names = []
    nt_lengths = []
    sequences = []
    metadata = []

    with open(csv_path) as f:
        reader = csv.reader(f, delimiter=",", quotechar='"')

        headers = next(reader)

        sequence_idx = None
        if has_sequences:
            sequence_idx = headers.index("sequence")
        names_idx = headers.index("name")

        headers.pop(names_idx)

        if has_sequences and names_idx < sequence_idx:
            sequence_idx -= 1

        headers.pop(sequence_idx)

        for row in reader:
            sequence_names.append(row.pop(names_idx))
            reference_genomes.append(reference_genome)

            if has_sequences:
                raw_seq = row.pop(sequence_idx)
                if "genbank://" in raw_seq:
                    raw_seq = raw_seq.split(":")

                    subset = None
                    if len(raw_seq) == 3:
                        subset = raw_seq[2]

                    ids = get_ids(raw_seq[1].replace("//", ""))
                    (nt_length, sequence) = fetch_sequence(ids, raw_seq, subset)
                    nt_lengths.append(str(nt_length))
                    sequences.append(sequence)
                else:
                    sequences.append(raw_seq)

            if len(metadata) == 0:
                metadata = [[] for _ in range(0, len(headers))]

            for (i, col) in enumerate(row):
                metadata[i].append(col)

    data.headers.extend(headers)
    data.columns = [reference_genomes, sequence_names, nt_lengths, sequences]
    data.columns.extend(metadata)
    return (data, config)
