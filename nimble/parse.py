import pysam

import pandas as pd

from Bio import SeqIO
from io import StringIO

from nimble.types import Data
from nimble.utils import get_library_name_from_filename


def parse_fasta(seq_path):

  data = Data()

  reference_name = get_library_name_from_filename(seq_path)

  for record in SeqIO.parse(seq_path, "fasta"):
    data.columns[0].append(reference_name)
    data.columns[1].append(record.id)
    data.columns[2].append(str(len(record)))
    data.columns[3].append(str(record.seq))

  return data


def parse_bam(seq_path):
  is_single_cell = False

  data = Data()

  # Add empty CellBarcode and UMI
  if is_single_cell:
    data.headers.extend(["UMI", "cell_barcode"])
    data.columns.extend([[], []])

  library_name = get_library_name_from_filename(seq_path)

  for read in pysam.AlignmentFile(seq_path, "rb"):
    seq = read.query_sequence

    data.columns[0].append(library_name)
    data.columns[1].append(read.reference_name)
    data.columns[2].append(len(seq))
    data.columns[3].append(seq)

  return data


def load_data_from_tsv(input_path):
  with open(input_path, "r") as f:
    metadata = [next(f)]

    str_rep = ""
    max_line_len = 0

    for line in f:
      csv_line = line.split("\t")[1].strip() + "," + line.split("\t")[0] + "\n"
      str_rep += csv_line + "\n"
      curr_line_len = len(csv_line.split(","))

      if(curr_line_len > max_line_len):
        max_line_len = curr_line_len

      metadata.append(line.split("\t")[1:])

  names = [i for i in range(0, max_line_len)]
  return (pd.read_csv(StringIO(str_rep), header=None, names=names), metadata)