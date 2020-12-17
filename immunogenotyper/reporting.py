import csv

def report(reference_path, results_path, output_path):
  with open(results_path, "r") as f:
    data = csv.reader(f, delimiter="\t")
    for row in data:
      print(row)