import pandas as pd
import numpy as np

from nimble.parse import parse_alignment_results
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


# Run the given method with the given value and return the data/metadata
def _filter(data, metadata, method, value):
  out_data = None

  if method == "minPct":
    out_data = _min_pct(data, value)
  elif method == "minCount":
    out_data = _min_count(data, value)
  else:
    raise ValueError("No such filter, " + method)

  return (out_data, metadata)


# API for this module. Can chain reports in an order provided by the methods list.
def report(methods, values, results_path, output_path):
  (data, metadata) = parse_alignment_results(results_path)

  for (method, value) in zip(methods, values):
    (data, metadata) = _filter(data, metadata, method, value)

  write_data_to_tsv(output_path, data, metadata)