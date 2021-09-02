import sys

USAGE_STRING = """
Documentation is at https://github.com/BimberLab/nimble/blob/main/README.md

usage:
  python -m nimble download <OPTIONAL:aligner_version>

  python -m nimble generate <fasta_path> | <csv_path> | <fasta_path> <csv_path> <output_path>

  python -m nimble align <reference.json> <input>...

  python -m nimble report <reference.json> <input.tsv> <output.tsv>

  python -m nimble filter <method> <method parameters...> <input.tsv> <output.tsv>
    Reporting methods:
      minPct <OPTIONAL:pct, default=0.01>
      minCount <OPTIONAL:count, default=5>

  python -m nimble help"""


def print_usage_and_exit():
    print(USAGE_STRING)
    sys.exit()
