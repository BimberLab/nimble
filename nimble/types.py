from enum import Enum


# Enum that specifies what type of data is in the output reference library
class DataType(str, Enum):
  FASTA = "fasta"
  BAM = "bam"
  SINGLECELL = "single-cell"
  

# Alignment config to be serialized and passed to the backend aligner
class Config():
  def __init__(self):
    self.score_threshold = 60
    self.score_filter = 25
    self.num_mismatches = 0
    self.discard_multiple_matches = False
    self.intersect_level = 0
    self.group_on = ""
    self.discard_multi_hits = 0
    self.require_valid_pair = False
    self.data_type = DataType.BAM
    self.filters = []


# Type to contain the actual sequence data/metadata. Can be arbitrarily modified at runtime to contain any metadata.
class Data():
  def __init__(self):
    self.headers = ["reference_genome", "sequence_name", "nt_length", "sequence"]
    self.columns = [[], [], [], []]