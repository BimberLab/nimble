import sys
import json
import argparse
from pathlib import Path
from Bio import SeqIO

USAGE_STRING = "usage: python main.py <fasta_path> <reference_json_path> <aligner_config_json_path>"

class Config():
  def __init__(self):
    self.score_threshold = 60
    self.percent_threshold = 0.05
    self.num_mismatches = 0
    self.discard_multiple_matches = False
    self.group_on = []


class Data():
  def __init__(self):
    self.headers = ["reference_genome", "nt_sequence", "nt_length"]
    self.rows = [[], [], []]


def main():
  seq_path = sys.argv[1]
  reference_output_path = sys.argv[2]
  config_output_path = sys.argv[3]

  reference_genome = Path(seq_path).stem.replace("_", " ")

  data = Data()

  for record in SeqIO.parse(seq_path, "fasta"):
    data.rows[0].append(reference_genome)
    data.rows[1].append(record.id,)
    data.rows[2].append(len(record))

  with open(reference_output_path, "w") as f:
    json.dump(data.__dict__, f, indent=2)

  with open(config_output_path, "w") as f:
    json.dump(Config().__dict__, f, indent=2)

  
if __name__ == "__main__":
  if len(sys.argv) != 4:
    print(USAGE_STRING)
    sys.exit()

  main()