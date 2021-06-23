from enum import Enum

class DataType(str, Enum):
  FASTA = "fasta"
  BAM = "bam"
  SINGLECELL = "single-cell"


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
    self.data_type = DataType.SINGLECELL


class Data():
  def __init__(self):
    self.headers = ["reference_genome", "sequence_name", "nt_length", "sequence"]
    self.columns = [[], [], [], []]
