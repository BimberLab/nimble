import pandas as pd
import pysam

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


# Given a path to a .bam, reads it and returns a (Data, Config) tuple with filled objects
def parse_bam(seq_path):
  data = Data()
  config = Config()
  config.data_type = DataType.BAM

  is_single_cell = False

  library_name = get_library_name_from_filename(seq_path)

  f = pysam.AlignmentFile(seq_path, "rb")

  for index, read in enumerate(f):
    cell_barcode = None
    molecular_barcode = None

    # Detect if any of these reads are 10x data. If they are, add the relevant columns and metadata
    try:
      cell_barcode = read.get_tag("CR")
      molecular_barcode = read.get_tag("UR")

      # Add empty CellBarcode and UMI
      if not is_single_cell:
        data.headers.extend(["cell_barcode", "molecular_barcode"])
        data.columns.extend([[], []])

        # Fill columns up to this read with None
        for _ in range(0, index):
          data.columns[4].append("null")
          data.columns[5].append("null")

        is_single_cell = True
        config.data_type = DataType.SINGLECELL
    except KeyError:
      pass

    seq = read.query_sequence

    data.columns[0].append(library_name)

    if read.reference_name is not None:
      data.columns[1].append(read.reference_name)
    else:
      data.columns[1].append("null")

    data.columns[2].append(str(len(seq)))
    data.columns[3].append(seq)

    if is_single_cell:
      if cell_barcode is not None:
        data.columns[4].append(cell_barcode)
      else:
        data.columns[4].append("null")

      if molecular_barcode is not None:
        data.columns[5].append(molecular_barcode)
      else:
        data.columns[4].append("null")


  return (data, config)


# Read data from the backend aligner's output format -- a TSV
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