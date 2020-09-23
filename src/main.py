import sys
import json
import argparse
from pathlib import Path
from Bio import SeqIO

USAGE_STRING = "usage: python main.py <fasta_path> <output_path>"

class Data():
  def __init__(self, reference_genome, nt_sequence, nt_length):
    self.reference_genome = reference_genome
    self.nt_sequence = nt_sequence
    self.nt_length = nt_length


def main():
  seq_path = sys.argv[1]
  output_path = sys.argv[2]

  reference_genome = Path(seq_path).stem.replace("_", " ")

  data = []

  for record in SeqIO.parse(seq_path, "fasta"):
    data.append(Data(reference_genome, record.id, len(record)).__dict__)

  with open(output_path, "w") as f:
    json.dump(data, f)
  
if __name__ == "__main__":
  if len(sys.argv) != 3:
    print(USAGE_STRING)
    sys.exit()

  main()