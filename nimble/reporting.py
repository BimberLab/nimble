import pandas as pd
import numpy as np

from nimble.parse import load_data_from_tsv
from nimble.utils import write_data_to_tsv


# Get the number of unique reference names in a readset
def _get_unique_references(data):
  return pd.unique(data[data.columns[1:]].values.ravel('K'))


# Filter reads via a threshold on the minimum percentage of total reads
def _min_pct(data, pct):
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


# Filter reads via a threshold on the minimum number of hits
def _min_count(data, count):
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


# API for this module
def report(method, value, results_path, output_path):
  (data, metadata) = load_data_from_tsv(results_path)

  out_data = None

  if method == "minPct":
    out_data = _min_pct(data, value)
  elif method == "minCount":
    out_data = _min_count(data, value)
  else:
    print("Incorrect format. Please see 'nimble help'.")

  write_data_to_tsv(output_path, out_data, metadata)