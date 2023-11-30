# nimble
nimble is a fast, accurate, and configurable RNA sequence alignment tool that executes lightweight alignments on arbitrary reference libraries. It uses pseudo-alignment to rapidly generate supplemental calls to complement a data pipeline's primary alignment. It does this with low overhead, making it possible to run supplemental alignments on almost any machine.


# Installation

nimble requires Python 3 and [samtools](http://www.htslib.org/). It supports the following operating systems:

- Ubuntu
- CentOS 7
- Manjaro

To install nimble, run the following command:

`pip install git+https://github.com/BimberLab/nimble`

There is also a [docker image](https://github.com/BimberLab/nimble/pkgs/container/nimble).


# Usage
Here's an example of generating a library and using it in an alignment, for both .bam and .fastq.gz inputs:
```
nimble generate --file lib.csv --opt-file lib.fasta --output_path lib.json

nimble align --reference lib.json,lib2.json,lib3.json --output out.tsv --input data.bam --alignment_path log.tsv.gz --log log.txt --num_cores 8 --strand_filter fiveprime

nimble align --reference lib.json,lib2.json,lib3.json --output out.tsv --input data_r1_fastq.gz data_r2_fastq.gz --alignment_path log.tsv.gz --log log.txt --num_cores 8 --strand_filter unstranded
```

A library is comprised of an input csv or fasta, or both.
`nimble align` can take one or more libraries via the `--reference` flag in a CSV string. Every other file will be generated per-library, and the library name will be appended to the filename root. For instance, in this case, the `--output` files will be `out-lib.tsv`, `out-lib2.tsv`, and `out-lib3.tsv`.

# Documentation

Detailed documentation can be found at the [wiki](https://github.com/BimberLab/nimble/wiki).

The source code for the backend aligner can be found [here](https://github.com/BimberLab/nimble-aligner).



# Issues

To report a bug, ask a question, or request support for a new operating system, please create an [issue](https://github.com/BimberLab/nimble/issues) in this repository.


# References

[Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.](http://bioinformatics.oxfordjournals.org/content/early/2015/10/06/bioinformatics.btv573.short?rss=1)

[rust-pseudoaligner](https://github.com/10XGenomics/rust-pseudoaligner)
