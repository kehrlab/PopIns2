# PopIns2

![GitHub release (latest by date)](https://img.shields.io/github/v/release/kehrlab/popins2)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/kehrlab/popins2/ci.yaml?branch=main)
![GitHub Repo stars](https://img.shields.io/github/stars/kehrlab/PopIns2?style=social)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4890793.svg)](https://doi.org/10.5281/zenodo.4890793)

Population-scale detection of non-reference sequence variants using colored de Bruijn Graphs

## Contents:
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Example](#example)
5. [Snakemake](#snakemake)
6. [Help](#help)
7. [Citation](#citation)

## Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 16.04 / 18.04, CentOS Linux 7.6 |
| C++14 capable compiler | g++ vers. 4.9.2, 5.5.0, 7.2.0 |
| [Bifrost](https://github.com/pmelsted/bfgraph) | vers. 1.0.4-ab43065 |
| [bwa](https://github.com/lh3/bwa) | vers. 0.7.15-r1140 |
| [samtools](https://github.com/samtools/samtools) | vers. 1.3, 1.5 |
| [sickle](https://github.com/najoshi/sickle) | vers. 1.33 |
| [gatb-minia-pipeline](https://github.com/Krannich479/gatb-minia-pipeline) | (*submodule; no need to install*) |
| [SeqAn](https://www.seqan.de/) | (*header library; no need to install*) |

Prior to the installation make sure your system meets all the requirements. For the default settings of PopIns2 a *Bifrost* installation with MAX_KMER_SIZE=64 is required. Presently, the conda package of Bifrost does not meet this requirement. If the executables of the software dependencies (bwa, samtools, sickle) are not accessible systemwide, you have to write the full paths to the executables into a configfile (see [Installation](#installation)). Submodules and header libraries come by default with the git clone, there is no need for a manual installation. For backward compatibility PopIns2 still offers to use the *Velvet assembler* (see [popins](https://github.com/bkehr/popins) for installation recommendation).
<br>
**Important update note**: With a release on April 28, 2022, the Bifrost API underwent major changes in the color implementation. The compatibility of PopIns2 with these later releases has not been tested thoroughly yet and might violate the objective of PopIns2. For the time being, using a Bifrost version prior to _commit 703be6d_ is recommended.

## Installation:

```
git clone --recursive https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

If the binaries of the software dependencies are not globally available on your system (e.g. by appending them to your `PATH`) you have to set the paths to the binaries within the *popins2.config* prior to executing `make`. After the compilation with `make` you should see the binary *popins2* in the main folder. The PopIns2 [Wiki](https://github.com/kehrlab/PopIns2/wiki) gathers known issues that might occur during installation or runtime.

## Usage:

PopIns2 is a program consisting of several submodules. The submodules are designed to be executed one after another and fit together into a consecutive workflow. To display the help page of a submodule type `popins2 <command> --help` as shown in the [help section](#help).

#### The assemble command
```
popins2 assemble [OPTIONS] sample.bam
```
The assemble command identifies reads without high-quality alignment to the reference genome, filters reads with poor base quality and assembles them into a set of contigs. The reads, given as BAM file, must be indexed by _bwa index_. Optionally, reads can be remapped to an additional reference FASTA before the filtering and assembly such that only the remaining reads without a high-quality alignment are further processed (e.g. useful for decontamination). The additional reference FASTQ must be indexed by _bwa index_ too.

#### The merge command
```
popins2 merge [OPTIONS] {-s|-r} DIR
```
\[Default\] The merge command builds a colored and compacted de Bruijn Graph (ccdbg) of all contigs of all samples in a given source directory _DIR_.
By default, the merge module finds all files of the pattern `<DIR>/*/assembly_final.contigs.fa`. To process the contigs of the [assemble command](#the-assemble-command) the __-r__ input parameter and graph simplification options __-d__ and __-i__ are highly recommended. Once the ccdbg is built, the merge module identifies paths in the graph and returns _supercontigs_.

```
popins2 merge [OPTIONS] -y GFA -z BFG_COLORS
```
An alternative way of providing input for the merge command is to directly pass a ccdbg. Here, the merge command expects a _GFA_ file and a _bfg_colors_ file, which is specific to the Bifrost. If you choose to run the merge command with a _pre_-built GFA graph, mind that you have to set the Algorithm options accordingly (in particular __-k__).

#### The contigmap command
```
popins2 contigmap [OPTIONS] SAMPLE_ID
```
The contigmap command maps all reads with low-quality alignments of a sample to the set of supercontigs using BWA-MEM. The mapping information is then merged with the reads' mates.

#### The place commands
```
popins2 place-refalign [OPTIONS]
popins2 place-splitalign [OPTIONS] SAMPLE_ID
popins2 place-finish [OPTIONS]
```
In brief, the place commands attempt to anker the supercontigs to the samples. At first, all potential anker locations from all samples are collected. Then prefixes/suffixes of the supercontigs are aligned to all collected locations. For successful alignments records are written to a VCF file. In the second step, all remaining locations are split-aligned per sample. Finally, all locations from all successful split-alignments are combined and added to the VCF file.

#### The genotype command
```
popins2 genotype [OPTIONS] SAMPLE_ID
```
The genotype command generates alleles (ALT) of the supercontigs with some flanking reference genome sequence. Then, the reads of a sample are aligned to ALT and the reference genome around the breakpoint (REF). The ratio of alignments to ALT and REF determines a genotype quality and a final genotype prediction per variant per sample.

## Example:

Test data for a minimum working example can be found at [zenodo](https://doi.org/10.5281/zenodo.4890793). A simple project structure for PopIns2 looks like

```
$ tree /path/to/your/project/
/path/to/your/project/
├── myFirstSample
│   ├── first_sample.bam
│   └── first_sample.bam.bai
├── mySecondSample
│   ├── second_sample.bam
│   └── second_sample.bam.bai
└── myThirdSample
    ├── third_sample.bam
    └── third_sample.bam.bai
```

and a simple workflow could look like

```
cd /path/to/your/project
ln -s /path/to/reference_genome.fa genome.fa
ln -s /path/to/reference_genome.fa.fai genome.fa.fai

popins2 assemble --sample sample1 /path/to/your/project/myFirstSample/first_sample.bam
popins2 assemble --sample sample2 /path/to/your/project/mySecondSample/second_sample.bam
popins2 assemble --sample sample3 /path/to/your/project/myThirdSample/third_sample.bam

popins2 merge -r /path/to/your/project -di

popins2 contigmap sample1
popins2 contigmap sample2
popins2 contigmap sample3

popins2 place-refalign
popins2 place-splitalign sample1
popins2 place-splitalign sample2
popins2 place-splitalign sample3
popins2 place-finish

popins2 genotype sample1
popins2 genotype sample2
popins2 genotype sample3
```

## Snakemake:

The workflow of PopIns2 can be effectively distributed among a HPC cluster environment. This [Github project](https://github.com/Krannich479/PopIns2_snakeproject) provides a template of a full PopIns2 workflow as individual cluster jobs using [Snakemake](https://snakemake.readthedocs.io/en/stable/), a Python-based workflow management tool.

## Help:

```
$ popins2 -h

Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs
================================================================

SYNOPSIS
    popins2 COMMAND [OPTIONS]

COMMAND
    assemble            Filter, clip and assemble unmapped reads from a sample.
    merge               Generate supercontigs from a colored compacted de Bruijn Graph.
    multik              Multi-k framework for a colored compacted de Bruijn Graph.
    contigmap           Map unmapped reads to (super-)contigs.
    place-refalign      Find position of (super-)contigs by aligning contig ends to the reference genome.
    place-splitalign    Find position of (super-)contigs by split-read alignment (per sample).
    place-finish        Combine position found by split-read alignment from all samples.
    genotype            Determine genotypes of all insertions in a sample.

VERSION
    0.12.0-a935f00, Date: on 2020-10-21 12:50:29

Try `popins2 COMMAND --help' for more information on each command.

```

For more troubleshooting, FAQs and tips about the usage of PopIns2 please have a look into the PopIns2 [Wiki](https://github.com/kehrlab/PopIns2/wiki).

## Citation:

Thomas Krannich, W Timothy J White, Sebastian Niehus, Guillaume Holley, Bjarni V Halldórsson, Birte Kehr. _Population-scale detection of non-reference sequence variants using colored de Bruijn graphs_. Bioinformatics 2022, 38(3):604–611. https://doi.org/10.1093/bioinformatics/btab749
