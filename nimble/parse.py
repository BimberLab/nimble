import json

import pandas as pd

from Bio import SeqIO
from io import StringIO

from nimble.types import Data, Config, DataType
from nimble.utils import get_library_name_from_filename


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