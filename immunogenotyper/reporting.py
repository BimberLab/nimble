from io import StringIO
import pandas as pd
from usage import print_usage_and_exit


def load_data(input_path):
  with open(input_path, "r") as f:
    next(f)

    str_rep = ""
    max_line_len = 0

    for line in f:
      csv_line = line.split("\t")[0]
      str_rep += csv_line + "\n"
      curr_line_len = len(csv_line.split(","))

      if(curr_line_len > max_line_len):
        max_line_len = curr_line_len

  names = [i for i in range(0, max_line_len)]
  return pd.read_csv(StringIO(str_rep), header=None, names=names)


def write_data(output_path, out_data):
  with open(output_path, "w") as f:
    pass


def min_pct(data, pct):
  if pct == None:
    pct = 0.01

  return ""


def min_count(data, value):
  if count == None:
    count = 5

  return ""


def min_pct_lineage(data, value):
  if pct == None:
    pct = 0.01

  return ""


def report(method, value, results_path, output_path):
  data = load_data(results_path)

  out_data = None

  if method == "minPct":
    out_data = min_pct(data, value)
  elif method == "minCount":
    out_data = min_count(data, value)
  elif method == "minPctLineage":
    out_data = min_pct_lineage(data, value)
  else:
    print_usage_and_exit()

  write_data(output_path, out_data)