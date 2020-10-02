import sys
import json
import argparse
from pathlib import Path
from Bio import SeqIO

USAGE_STRING = """usage:
  python main.py <fasta_path> <reference_json_path> <aligner_config_json_path>

  OR

  python main.py compile <reference_json_path> <aligner_config_json_path> <compiled_json_path>
  """

class Config():
  def __init__(self):
    self.score_threshold = 60
    self.score_filter = 25
    self.num_mismatches = 0
    self.discard_multiple_matches = False
    self.group_on = ""


class Data():
  def __init__(self):
    self.headers = ["reference_genome", "nt_sequence", "nt_length"]
    self.columns = [[], [], []]


# Take human-editable config json files and compile into a single minified file for the aligner
def compile_config():
  reference_output_path = sys.argv[2]
  config_output_path = sys.argv[3]
  compiled_json_path = sys.argv[4]

  with open(reference_output_path, "r") as ref, open(config_output_path, "r") as conf, open(compiled_json_path, "w") as comp:
    reference = json.load(ref)
    config = json.load(conf)
    json.dump([config, reference], comp)


# Generate and write human-editable config json files to disk
def create_editable_config():
  seq_path = sys.argv[1]
  reference_output_path = sys.argv[2]
  config_output_path = sys.argv[3]

  # Generate reference genome name by getting the filename and prettifying it
  reference_genome = Path(seq_path).stem.replace("_", " ")

  data = Data()

  for record in SeqIO.parse(seq_path, "fasta"):
    data.columns[0].append(reference_genome)
    data.columns[1].append(record.id,)
    data.columns[2].append(str(len(record)))

  with open(reference_output_path, "w") as f:
    json.dump(data.__dict__, f, indent=2)

  with open(config_output_path, "w") as f:
    json.dump(Config().__dict__, f, indent=2)


 
if __name__ == "__main__":
  print_usage_and_exit = False

  if len(sys.argv) == 1:
      print_usage_and_exit = True
  elif sys.argv[1] == "compile":
    if len(sys.argv) != 5:
      print_usage_and_exit = True
    else:
      compile_config()
  else:
    if len(sys.argv) != 4:
      print_usage_and_exit = True
    else:
      create_editable_config()

  if print_usage_and_exit:
    print(USAGE_STRING)
    sys.exit()
