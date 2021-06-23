import pandas as pd
import numpy as np

from nimble.parse import load_data_from_tsv
from nimble.utils import write_data_to_tsv


def _get_unique_references(data):
  return pd.unique(data[data.columns[1:]].values.ravel('K'))


def min_pct(data, pct):
  if pct == None:
    pct = 0.01

  num_reads_total = data[0].sum()
  references = _get_unique_references(data)
  references_to_drop = []

  for reference in references:
    num_reads = data[data.apply(lambda row: reference in row.values, axis=1)][0].sum()

    if (num_reads / num_reads_total < pct):
      references_to_drop.append(reference)

  for reference in references_to_drop:
    data = data.replace(reference, np.nan)

  return data


def min_count(data, count):
  if count == None:
    count = 5

  references = _get_unique_references(data)
  references_to_drop = []

  for reference in references:
    num_reads = data[data.apply(lambda row: reference in row.values, axis=1)][0].sum()

    if (num_reads < count):
      references_to_drop.append(reference)

  for reference in references_to_drop:
    data = data.replace(reference, np.nan)

  return data


def min_pct_lineage(data, value):
  if pct == None:
    pct = 0.01

  return ""


def report(method, value, results_path, output_path):
  (data, metadata) = load_data_from_tsv(results_path)

  out_data = None

  if method == "minPct":
    out_data = min_pct(data, value)
  elif method == "minCount":
    out_data = min_count(data, value)
  elif method == "minPctLineage":
    out_data = min_pct_lineage(data, value)
  else:
    print("Incorrect format. Please see 'nimble help'.")

  write_data_to_tsv(output_path, out_data, metadata)