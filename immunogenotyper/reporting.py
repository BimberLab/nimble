import csv


def simple_filter(data):
  pass


def report(method, results_path, output_path):
  with open(results_path, "r") as f:
    data = csv.reader(f, delimiter="\t")

    if method == "simple":
      simple_filter(data)
    else:
      print("\nNo such reporting method. Please choose from the following list (documentation is at https://github.com/devsebb/ImmunoGenotyperLibraryifier/blob/main/README.md):")
      print("\tsimple")
      print("")