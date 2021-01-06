import sys

USAGE_STRING = """
Documentation is at https://github.com/devsebb/ImmunoGenotyperLibraryifier/blob/main/README.md

usage:
  immunogenotyper download <OPTIONAL:aligner_version>

  immunogenotyper generate <fasta_path> <reference_json_path> <aligner_config_json_path>

  immunogenotyper compile <reference_json_path> <aligner_config_json_path> <compiled_json_path>

  immunogenotyper align <reference.json> <reference.fasta> <input>...

  immunogenotyper report <method> <method parameters...> <input.tsv> <output.tsv>
    Reporting methods:
      minPct <OPTIONAL:pct, default=0.01>
      minCount <OPTIONAL:count, default=5>
      minPctLineage <OPTIONAL:pct, default=0.01>

  immunogenotyper help"""

def print_usage_and_exit():
  print(USAGE_STRING)
  sys.exit()