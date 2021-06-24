import sys
import distro
import pandas as pd

from pathlib import Path


# Generate reference genome name by getting the filename and prettifying it
def get_library_name_from_filename(seq_path):
  return Path(seq_path).stem.replace("_", " ")


# Get the name of a Github release given the target platform
def get_exec_name_from_platform():
  exec_name = ""

  if sys.platform == "win32":
    exec_name = "windows.exe"
  elif sys.platform == "linux":
    # Detect Centos7 vs other distros
    if distro.id() == "centos":
      exec_name = "CentOS.out"
    elif distro.id() == "manjaro":
      exec_name = "Manjaro.out"
    else:
      exec_name = "Ubuntu.out"
  elif sys.platform == "darwin":
    exec_name = "Mac.app"
  else:
    print("Error -- platform " + sys.platform + " does not have a compatible release on GitHub.")
    sys.exit()

  return exec_name


# Serialize results of an alignment/filter to output TSV
def write_data_to_tsv(output_path, out_data, metadata):
  out_data = out_data.drop(0, axis=1)

  with open(output_path, "w") as f:
    metadata = iter(metadata)
    f.write(next(metadata))

    for row in out_data.apply(lambda row: row.values[~pd.isna(row.values)], axis=1):
      line_metadata = "\t".join([elem for elem in next(metadata)])
      
      if len(row) > 0:
        f.write(",".join([reference for reference in row]) + "\t" + line_metadata)