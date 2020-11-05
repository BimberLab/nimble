# Immunogenotyper
Immunogenotyper is a developing program for aligning Bulk-seq and 10x RNA data against arbitrary reference genomes. This repo contains the end-user command line interface and related documentation -- for details about the backend aligner, go [here](https://github.com/devsebb/ImmunoGenotyper).


# Installation
To install Immunogenotyper, ensure you have a version of Python >= 3.6. Then, run the following command:

`pip install git+https://github.com/devsebb/ImmunoGenotyperLibraryifier`



# Usage

To perform an alignment of RNA data against your reference library, do the following:

1. Run `python -m immunogenotyper generate <reference-library.fasta> <output-reference.json> <output-config.json>`. `<reference-library.fasta>` should be the file containing your reference library sequences, `<output-reference.json>` is the name/path of the file you want to create containing the reference library metadata, and `<output-config.json>` is the name/path of the file you want to create containing a default aligner configuration.

2. Edit the `<output-reference.json>` and `<output-config.json>` files created by the previous step. For instance, if you want to add allele-level metadata like lineage or locus to the reference library, add them as a column to the `<output-reference.json>` now. The `<output-config.json>` contains several values for configuring the aligner, so edit those if necessary as well.

3. Run `python -m immunogenotyper compile <output-reference.json> <output-config.json> <output-compiled.json>`. This takes the two .json files created by the previous steps and produces a combined file at the location you specify with `<output-compiled.json>`.

4. Run `python -m immunogenotyper align <output-compiled.json> <reference-library.fasta> <input-r1.fastq> <input-r2.fastq:OPTIONAL>`. `<output-compiled.json>` is the file created by the `compile` command in the previous step, and `<reference-library.fasta>` is your original reference library sequence file. `<input-r1.fastq>` and `<input-r2.fastq>` are paths to your input bulk-seq data (`<input-r2.fastq>` is optional). This command will produce a file, `results.tsv`, with the results of the alignment.


The `align` command takes an optional flag, `--debug-reference <REFERENCE_NAME>` which produces a `debug.tsv` file containing the sequences and scores of each read that matched the reference `<REFERENCE_NAME>`.

The `download` command downloads the most recent aligner [release](https://github.com/devsebb/ImmunoGenotyper/releases). You can optionally specify a release version like so:
`python -m immunogenotyper download <version>` where `<version>` is a release tag like "v0.0.1" or "v0.0.1-beta.1". Note that this command downloads a file called `aligner` to the directory you're running the command in -- be careful not to overwrite an existing file.

To get help information, run `python -m immunogenotyper` or `python -m immunogenotyper help`.
