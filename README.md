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



# JSON Format
The `generate` command creates two .json files. One file contains the reference metadata, and the other contains the configuration for the aligner.


## Reference Metadata
The reference metadata file follows this format:
```
{
  "headers": ["reference_genome", "nt_sequence", "nt_length", ...]
  "columns": [[...], [...], [...]]
}
```

This file contains a `headers` field and a `columns` field. `headers` is an array of strings that corrospond to the matching column in the `columns` field. The aligner must have at least a `reference_genome` header, an `nt_sequence` header, and an `nt_length` header.

The `columns` field is a multidimentional array of strings. Each sub-array corrosponds to a header in the `headers` field.

To add another header/column pair (e.g. to add per-allele lineage or locus information), add a string to the `headers` array and add a column to the corrosponding index in the `columns` field.


## Aligner Configuration
The aligner configuration file follows this format:

```
{
  "score_threshold": number,
  "score_filter": number,
  "num_mismatches": number,
  "discard_multiple_matches": boolean,
  "intersect_level: number",
  "group_on": string
}
```

* `score_threshold`: controls the score an alignment needs to reach to be considered a match. For perfect matches, set this value equal to the length of the reads being aligned to the reference library.

* `score_filter`: sets a lower boundary on the number of matches needed on a reference before it is reported. For instance, if you set `"score_filter": 25`, no reference with less than 25 matches will be reported in the output.

* `num_mismatches`: sets the allowable number of mismatches during alignment.

* `discard_multiple_matches`: flag for whether a read that matches multiple references should be counted. If `true`, a read that matches multiple references will count toward the scores of all of those references. If `false`, the read's matches are discarded.

* `intersect_level`: controls logic behind how to count matches during alignment. There are three intersect levels. `intersect_level: 0` takes the best matches from either the read or reverse read, determined by alignment score. `intersect_level: 1` takes the intersection between the read and reverse read -- if there is no intersection, it defaults to the best match. `intersect_level: 2` takes the intersection and reports no match if there is no intersection.

* `group_on`: if this is set to the name of a header in the reference metadata file, the output `results.tsv` will be filtered to that level of specificity. For instance, if you've added a column with lineage information under a header called "lineage", sestting `"group_on": "lineage"` will report lineage-level information, rather than the default case of allele-level information. If a single read matches onto the `group_on` category more than once during alignment (for instance, if a read matches multiple alleles in the same lineage and you're grouping on lineage), it will only count as one match. If `group_on` is unset, allele-level information is returned.
