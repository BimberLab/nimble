import csv
from usage import print_usage_and_exit


def minPct(data, pct):
  if pct == None:
    pct = 0.01


def minCount(data, value):
  if count == None:
    count = 5


def minPctLineage(data, value):
  if pct == None:
    pct = 0.01


def report(method, value, results_path, output_path):
  with open(results_path, "r") as f:
    data = csv.reader(f, delimiter="\t")

    if method == "minPct":
      minPct(data, value)
    elif method == "minCount":
      minCount(data, value)
    elif method == "minPctLineage":
      minPctLineage(data, value)
    else:
      print_usage_and_exit()