# Welcome to MkDocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

## Commands

Assemble a genome!

Learn how to create and assess genome assemblies using the powerful combination of Nanopore and Illumina reads

This tutorial explores how long and short read data can be combined to produce a high-quality ‘finished’ bacterial genome sequence. Termed ‘hybrid assembly’, we will use read data produced from two different sequencing platforms, Illumina (short read) and Oxford Nanopore Technologies (long read), to reconstruct a bacterial genome sequence.

In this tutorial we will perform ‘de novo assembly’. De novo assembly is the process of assembling a genome from scratch using only the sequenced reads as input - no reference genome is used. This approach is common practise when working with microorganisms, and has seen increasing use for eukaryotes (including humans) in recent times.

Long reads can be used together with short reads to produce a high-quality assembly. Nanopore long reads (commonly >40,000 bases) can fully span repeats, and reveal how all the genome fragments should be arranged. Therefore, while long reads will provide the general structure of the genome, short reads will provide that high base-level accuracy needed to close a genome.

**Data:** Nanopore reads, Illlumina reads, bacterial organism (_Vibrio parahaemolyticus_) reference genome

**Tools:** Flye, Pilon, Unicycler, Quast, BUSCO, BWA, SAMTOOLS, FASTQC, Fastp,


## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
