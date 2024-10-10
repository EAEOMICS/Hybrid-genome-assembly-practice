# Welcome to MkDocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

## Commands

Assemble a genome!

Learn how to create and assess genome assemblies using the powerful combination of Nanopore and Illumina reads

This tutorial explores how long and short read data can be combined to produce a high-quality ‘finished’ bacterial genome sequence. Termed ‘hybrid assembly’, we will use read data produced from two different sequencing platforms, Illumina (short read) and Oxford Nanopore Technologies (long read), to reconstruct a bacterial genome sequence.

In this tutorial we will perform ‘de novo assembly’. De novo assembly is the process of assembling a genome from scratch using only the sequenced reads as input - no reference genome is used. This approach is common practise when working with microorganisms, and has seen increasing use for eukaryotes (including humans) in recent times.

Long reads can be used together with short reads to produce a high-quality assembly. Nanopore long reads (commonly >40,000 bases) can fully span repeats, and reveal how all the genome fragments should be arranged. Therefore, while long reads will provide the general structure of the genome, short reads will provide that high base-level accuracy needed to close a genome.

**Data:** Nanopore reads, Illlumina reads, bacterial organism (_Vibrio parahaemolyticus_) reference genome

**Tools:** `Flye`, `Pilon`, `Unicycler`, `Quast`, `Busco`, `BWA`, `Samtools`,`FastQC`, `Fastp`,`nanoplot`,`chopper`,`porechop`,`cutadapt`,`trimmomatic`

## Section 1: Read inspection and QC

In this section we will import and perform quality control (QC) on our data. 

Today we will use 4 pieces of data - **2 short read sets, 1 long read set, and a reference genome** to compare our assembly with. 

<br>

### Getting the data

To get the data we will import our reads directly from the **NCBI SRA database** (Sequence Read Archive) which is the largest publicly available repository of high throughput sequencing data

To download the data we will use the `fasterq-dump` tool from the SRA-toolkit which should be already installed in your conda environment.
We know that the accession ID of our samples is **SRR10022815** and **SRR10022816**. 

Before we start typing and running commands, is important to generate a new directory where all the files are going to be saved and organized.
In our case we are going to create a directory called `hybrid_assembly` and a subdirectory called `data` in this subdirectory is where we are going to save our reads. **Do not move them at any moment**

```bash
mkdir -p hybrid_assembly/data
cd hybrid_assembly/data
fasterq-dump {accesion_ID} #repeat for both accesion IDs
```

### Read inspection

Often, it is prudent to first assess the quality of our read sets. For the short reads, we are concerned with base quality, sequence duplication, and presence of adapter sequences. For nanopore, we want to know about the length and quality distribution of reads, as these may both be highly variable. 

`FastQC` creates summary reports for short read data. We will use this tool twice - once for each Illumina read set. We can then use a tool called `MultiQC` to combine these reports for easy viewing. 

For Nanopore data, `NanoPlot` is a great option. It creates plots which aim to summarise the length and quality distribution of long read sets. 

Depending on these summaries, we may choose to perform a QC step to remove any poor quality reads before proceeding. 


**Run FastQC**

As said before, we will use `FastQC` to see the quality of our **Illumina** reads. Once again, remember **Organization is key**, therefore we are going to create a new directory where we are going to save all our Quality Control outputs
```{Bash}
mkdir -p qc/illumina_raw
cd qc/illumina_raw

fastqc ../../data/SRR10022816_1.fastq ../../data/SRR10022816_2.fastq -o ./
```

FastQC produces two outputs - 'RawData', and 'Webpage'. Typically, the webpage is for human viewing, and the RawData can be given to other programs, such as MultiQC.

Let's see how the QC went on the Illumina Reads
<br>
<img src="assets/illuminaR1_raw.png">
<img src="assets/illuminaR2_raw.png">
<br>

<details>
<summary>Question 1(click to reveal)</summary>
What do you think about them? Do you think they have enough quality? Let's discuss, take your time to inspect the whole html
</details>


**Run Fastp**


**Run NanoPlot**
* Find MultiQC in the tools panel. It is listed as '**NanoPlot** Plotting suite for Oxford Nanopore sequencing data and alignments'
* Parameters:
    * ***Type of the file(s) to work on***
        * ***files*** nanopore_reads.fastq
* Leave all else default and execute the program. 

<br>

NanoPlot produces 5 outputs, but we are only interested in the 'HTML report' output. View this file by clicking the eye icon on this history item.

The NanoPlot HTML report includes a table, followed by a number of plots. The table provides a summary of the read set. The main statistics we will look at are ***Median read length***, ***Median read quality***, and ***Number of reads***. 

<br>

<img src="media/nanoplot_table.PNG" width="400">

<br>
