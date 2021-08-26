import csv
from Bio import Entrez

from nimble.types import Data, Config, DataType
from nimble.utils import get_library_name_from_filename

def get_remote_lib(tsv_path):
    Entrez.email = "benjamse@ohsu.edu"

    data = Data()
    config = Config()
    config.data_type = DataType.FASTA

    terms = []
    ids = []

    reference_genome = get_library_name_from_filename(tsv_path)
    reference_genomes = []
    sequence_names = []
    nt_lengths = []
    sequences = []

    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')

        next(reader)

        for row in reader:
            sequence_names.append(row[0])
            terms.append(row[1])

    for term in terms:
        ids += get_ids(term)

    handle = Entrez.efetch(db="nucleotide", id=ids, retmode="text", rettype="fasta")
    record = handle.read()

    lines = iter(record.splitlines())
    next(lines)

    curr_seq = ""

    for line in lines:
        if line and line[0] == ">":
            reference_genomes.append(reference_genome)
            nt_lengths.append(str(len(curr_seq)))
            sequences.append(curr_seq)
            curr_seq = ""
        else:
            curr_seq += line

    nt_lengths.append(str(len(curr_seq)))
    sequences.append(curr_seq)

    data.columns = [reference_genomes, sequence_names, nt_lengths, sequences]
    return (data, config)


def get_ids(term):
    ids = []
    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    ids += record["IdList"]
    return ids
