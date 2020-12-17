import sys
from sys import platform
import requests 
import os
import subprocess
import stat
import json
import argparse
from pathlib import Path
from Bio import SeqIO

from reporting import report

USAGE_STRING = """usage:
  immunogenotyper download <OPTIONAL:aligner_version>

  immunogenotyper generate <fasta_path> <reference_json_path> <aligner_config_json_path>

  immunogenotyper compile <reference_json_path> <aligner_config_json_path> <compiled_json_path>

  immunogenotyper align <reference.json> <reference.fasta> <input>...

  immunogenotyper report <reference.json> <results.tsv> <output.tsv>

  immunogenotyper help
  """

class Config():
  def __init__(self):
    self.score_threshold = 60
    self.score_filter = 25
    self.num_mismatches = 0
    self.discard_multiple_matches = False
    self.intersect_level = 0
    self.group_on = ""


class Data():
  def __init__(self):
    self.headers = ["reference_genome", "nt_sequence", "nt_length"]
    self.columns = [[], [], []]


# Take human-editable config json files and compile into a single minified file for the aligner
def compile_config(reference_output_path, config_output_path, compiled_json_path):
  with open(reference_output_path, "r") as ref, open(config_output_path, "r") as conf, open(compiled_json_path, "w") as comp:
    reference = json.load(ref)
    config = json.load(conf)
    json.dump([config, reference], comp)


# Generate and write human-editable config json files to disk
def create_editable_config(seq_path, reference_output_path, config_output_path):
  # Generate reference genome name by getting the filename and prettifying it
  reference_genome = Path(seq_path).stem.replace("_", " ")

  data = Data()

  # Fill data structure containing the reference library
  for record in SeqIO.parse(seq_path, "fasta"):
    data.columns[0].append(reference_genome)
    data.columns[1].append(record.id,)
    data.columns[2].append(str(len(record)))

  # Write reference and default config to disk
  with open(reference_output_path, "w") as f:
    json.dump(data.__dict__, f, indent=2)

  with open(config_output_path, "w") as f:
    json.dump(Config().__dict__, f, indent=2)


# Get the name of a Github release given the target platform
def get_exec_name_from_platform():
  exec_name = ""

  if sys.platform == "win32":
    exec_name = "windows.exe"
  elif sys.platform == "linux":
    exec_name = "Linux.x86"
  elif sys.platform == "darwin":
    exec_name = "Mac.app"
  else:
    print("Error -- platform " + sys.platform + " does not have a compatible release on GitHub.")
    sys.exit()

  return exec_name


# Given the name of a release, download the platform-specific executable from that release.
# Given no name, default to the most recent release.
def download_aligner(release):
  exec_name = get_exec_name_from_platform()

  url = ""
  if len(release) == 1:
    url = "https://github.com/devsebb/ImmunoGenotyper/releases/download/" + release[0] + "/" + exec_name
  else:
    url = "https://github.com/devsebb/ImmunoGenotyper/releases/latest/download/" + exec_name

  # Download aligner
  r = requests.get(url)

  # If we successfully downloaded it, write the file to disk and give it execution permissions if necessary
  if r.status_code == 200:
    with open('aligner', 'wb') as f:
      f.write(r.content)
    if exec_name == "Linux.x86" or exec_name == "Mac.app":
      st = os.stat('aligner')
      os.chmod('aligner', st.st_mode | stat.S_IEXEC)
  else:
    print("Error -- could not download aligner, status code " + str(r.status_code))
    sys.exit()


# Check if the aligner exists -- if it does, call it with the given parameters.
def align(param_list):
  path = os.path.abspath("aligner")
  
  if os.path.exists(path):
    subprocess.call([path] + param_list)
  else:
    print("Error -- no aligner found. Please run the 'immunogenotyper download' command.\n")
    sys.exit()


if __name__ == "__main__":
  if len(sys.argv) == 1: # Ensure we can index sys.argv[1]
    print(USAGE_STRING)
    sys.exit()
  elif sys.argv[1] == "download" and len(sys.argv) <= 3:
    download_aligner(sys.argv[2:])
  elif sys.argv[1] == "generate" and len(sys.argv) == 5:
    create_editable_config(sys.argv[2], sys.argv[3], sys.argv[4])
  elif sys.argv[1] == "compile" and len(sys.argv) == 5:
    compile_config(sys.argv[2], sys.argv[3], sys.argv[4])
  elif sys.argv[1] == "align":
    align(sys.argv[2:])
  elif sys.argv[1] == "report" and len(sys.argv) == 5:
    report(sys.argv[2], sys.argv[3], sys.argv[4])
  else:
    print(USAGE_STRING)
    sys.exit()
