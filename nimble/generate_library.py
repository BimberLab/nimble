import json
from nimble.parse import parse_fasta, parse_csv

# Generate and write nimble library json files to disk. Input data is a CSV, FASTA, or both
def generate(file, opt_file, output_path):

    (data, config, is_csv_req) = process_file(file, opt_file)
    (data_opt, config_opt, is_csv_opt) = process_file(opt_file, file)

    final_config = config
    if data_opt != None and is_csv_opt:
        final_config = config_opt

    final_data = None
    if data_opt != None:
        if is_csv_req:
            final_data = collate_data(data_opt, data)
        elif is_csv_opt:
            final_data = collate_data(data, data_opt)
    else:
        final_data = data

    print("Filtered " + str(low_complexity_filter_amount) + " base pairs from reference library.")

    # Write reference and default config to disk
    with open(output_path, "w") as f:
        json.dump([ final_config.__dict__, final_data.__dict__], f, indent=2)


# Parse a lone FASTA/lone CSV/CSV-FASTA tuple. If there's a lone FASTA, generate a simple library.
# If there's a lone CSV, assume it has sequence information or a genbank link.
# If there is not a lone CSV, assume it contains the metadata and that the sequence information is contained in the FASTA.
def process_file(file, paired_file):
    data = None
    config = None
    is_csv = False

    if file:
        if pathlib.Path(file).suffix == ".fasta":
            (data, config) = parse_fasta(file)
        elif pathlib.Path(file).suffix == ".csv" and paired_file:
            (data, config) = parse_csv(file, False)
            is_csv = True
        elif pathlib.Path(file).suffix == ".csv" and not paired_file:
            (data, config) = parse_csv(file, True)
            is_csv = True

    return (data, config, is_csv)