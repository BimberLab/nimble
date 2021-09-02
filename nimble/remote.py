import os
from Bio import Entrez

from nimble.types import Data, Config, DataType
from nimble.utils import get_library_name_from_filename


Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")


def fetch_sequence(ids, string_id, subset):
    if len(ids) > 1:
        print("Error -- attempt to fetch sequence with multiple ids: " + string_id)
        sys.exit()

    if len(ids) == 0:
        print("Error -- attempt to fetch sequence with no ids: " + string_id)
        sys.exit()

    handle = Entrez.efetch(db="nucleotide", id=ids, retmode="text", rettype="fasta")
    record = handle.read()

    lines = iter(record.splitlines())
    next(lines)

    seq = ""

    for line in lines:
        if line and line[0] == ">":
            reference_genomes.append(reference_genome)
            nt_lengths.append(str(len(seq)))
            sequences.append(seq)
            seq = ""
        else:
            seq += line

    if subset:
        subset = subset.split("-")
        seq = seq[int(subset[0]):int(subset[1])]
    return (len(seq), seq)


def get_ids(term):
    ids = []
    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    ids += record["IdList"]
    return ids
