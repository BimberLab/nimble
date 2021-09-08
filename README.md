# nimble
nimble is a fast, accurate, and configurable RNA sequence aligner that executes lightweight alignments on arbitrary reference libraries. It uses pseudo-alignment to rapidly generate supplemental calls to complement a data pipeline's primary alignment. It does this with low overhead, making it possible to run supplemental alignments on almost any machine.


# Installation

nimble requires Python 3 and [samtools](http://www.htslib.org/). It supports the following operating systems:

- Windows
- MacOS
- Ubuntu
- CentOS 7
- Manjaro

To install nimble, run the following command:

`pip install git+https://github.com/BimberLab/nimble`

There is also a [docker image](https://github.com/BimberLab/nimble/pkgs/container/nimble).


# Usage

For usage documentation, refer to the [quickstart guide](https://github.com/BimberLab/nimble/wiki/Quickstart).

There are also [usage examples](https://github.com/BimberLab/nimble/wiki/Example-Data-Analysis).


# Documentation

Detailed documentation can be found at the [wiki](https://github.com/BimberLab/nimble/wiki).

The source code for the backend aligner can be found [here](https://github.com/BimberLab/nimble-aligner).



# Issues

To report a bug, ask a question, or request support for a new operating system, please create an [issue](https://github.com/BimberLab/nimble/issues) in this repository.


# References

[Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.](http://bioinformatics.oxfordjournals.org/content/early/2015/10/06/bioinformatics.btv573.short?rss=1)

[rust-pseudoaligner](https://github.com/10XGenomics/rust-pseudoaligner)
