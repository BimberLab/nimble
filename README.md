# nimble
nimble is a lightweight tool designed to supplement standard RNA-seq and scRNA-seq pipelines by providing accurate gene quantification with customizable reference libraries. Traditional transcript quantification methods align reads to a single reference genome and apply a uniform feature-calling approach across all genes. While effective in most cases, this method can introduce systematic inaccuracies, particularly for complex or highly variable gene families such as MHC and KIR, where standard approaches may result in misalignment, lost counts, or incorrect feature assignments. Nimble addresses these gaps by allowing you to define custom gene spaces with configurable scoring criteria, ensuring more accurate and biologically relevant quantification.

![nimble_scrnaseq_pipeline](https://github.com/user-attachments/assets/baaa4015-0f4b-4c6e-bc9f-8155af4cd72b)

# Documentation
For installation instructions, [see below](#installation). Detailed documentation, including [usage instructions](https://github.com/BimberLab/nimble/wiki/Quickstart) and [a custom HLA reference vignette](https://github.com/BimberLab/nimble/wiki/HLA-Vignette), can be found at the [wiki](https://github.com/BimberLab/nimble/wiki).

# Support
Please [open an issue](https://github.com/BimberLab/nimble/issues/new) on this repository for feature requests, support for additional operating systems, bugfixes, or questions.

# Installation
The best way to install nimble is via our [Docker image](https://github.com/BimberLab/nimble/pkgs/container/nimble).

To install the latest nimble Docker image, run:

`docker pull ghcr.io/bimberlab/nimble:latest`

Alternatively, you can install nimble through pip. Nimble requires Python 3 and [samtools](http://www.htslib.org/). Currently, we support the following operating systems:

- MacOS
- Linux distributions with [musl](https://musl.libc.org/) support, i.e. Alpine and Debian/Ubuntu with the `musl` package, etc.

To install nimble with pip, run:

`pip install git+https://github.com/BimberLab/nimble`

Once you have nimble installed, proceed to our [usage documentation](https://github.com/BimberLab/nimble/wiki/Quickstart).

# Other resources
- The [nimbler](https://github.com/BimberLab/nimbleR) package contains R routines for integrating nimble data into a Seurat object.

- The source code for the backend aligner can be found [here](https://github.com/BimberLab/nimble-aligner).

# References and Credits
[rust-pseudoaligner](https://github.com/10XGenomics/rust-pseudoaligner)

[Samtools and HTSlib](www.htslib.org)

[Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.](http://bioinformatics.oxfordjournals.org/content/early/2015/10/06/bioinformatics.btv573.short?rss=1)

[Cock, P.J. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), pp.1422–1423.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2682512/)
