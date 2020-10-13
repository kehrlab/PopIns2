# PopIns2
Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

## Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 16.04 / 18.04, CentOS Linux 7.6 |
| C++14 capable compiler | g++ vers. 4.9.2, 5.5.0, 7.2.0 |
| [SeqAn](https://www.seqan.de/) | vers. 2.2.0 |
| [Bifrost](https://github.com/pmelsted/bfgraph) | vers. 1.0.4-ab43065 |
| [bwa](https://github.com/lh3/bwa) | vers. 0.7.15-r1140 |
| [samtools](https://github.com/samtools/samtools) | vers. 1.3, 1.5 |
| [sickle](https://github.com/najoshi/sickle) | vers. 1.33 |
| [gatb-minia-pipeline](https://github.com/Krannich479/gatb-minia-pipeline) | (*submodule; no need to install*) |

Prior to the installation make sure your system meets all the requirements. For the default settings of PopIns2 a *Bifrost* installation with MAX_KMER_SIZE=64 is required. If the executables of the software dependencies (bwa, samtools, sickle) are not accessible systemwide, you have to write the full paths to the executables into a configfile (see [Installation](#installation)). Submodules come by default with the git clone, there is no need for a manual installation. For backward compatibility PopIns2 still offers to use the *Velvet assembler* (see [popins](https://github.com/bkehr/popins) for installation recommendation).

## Installation:

```
git clone --recursive https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

If the binaries of the software dependencies are not globally available on your system (e.g. by appending them to your `PATH`) you have to set the paths to the binaries within the *popins2.config* prior to executing `make`. After the compilation with `make` you should see the binary *popins2* in the main folder.

## Usage:

PopIns2 is a program consisting of several submodules. To display the help page of a submodule type `popins2 <command> --help` as shown in the [help section](#help). 

#### The assemble command
```
popins2 assemble [OPTIONS] sample.bam
```
The assemble command identifies reads without high-quality alignment to the refence genome, filters reads with poor base quality and assembles them into a set of contigs. Optionally, reads can be remapped to an additional reference FASTA before the filtering and assembly such that only the remaining reads without a high-quality alignment are further processed (e.g. useful for decontamination). The additional reference FASTQ must be indexed using _bwa_. The assembly of contigs can be turned off.

#### The merge command
```
popins2 merge [OPTIONS] {-s|-r} DIR
```
The merge command builds a colored and compacted de Bruijn Graph (ccdbg) of all contigs of all samples in a given source directory _DIR_. 
By default, the merge module finds all files of the pattern `<DIR>/*/assembly_final.contigs.fa`. To process the contigs of the [assemble command](#the-assemble-command) the __-r__ input parameter is recommended. Once the ccdbg is built, the merge module identifies paths in the graph and returns _supercontigs_. \
 \
```
popins2 merge [OPTIONS] -y GFA -z BFG_COLORS
```
An alternative way of providing input for the merge command is to directly pass a ccdbg, e.g. as returned from the [multik command](#the-multik-command).
Here, the merge command expects a _GFA_ file and a _bfg_colors_ file, which is specific to the Bifrost library. If you choose to run the merge command with a _pre_-built GFA graph, mind that you have to set the Algorithm options accordingly (in particular __-k__).

#### The multik command
```
TODO
```
TODO

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
    0.11.0-0a8d447, Date: on 2020-10-08 17:06:03

Try `popins2 COMMAND --help' for more information on each command.

```

## Troubleshooting / FAQ:

- **Q1:** Where do I install _SeqAn_ and _Bifrost_ for PopIns2? **A:** The C++ libraries SeqAn and Bifrost have to be found by your compiler, i.e. they should usually be located in your system's `local/include` and `local/lib` folders. The respective website and github page (see [Requirements](#requirements)) provide more details.

- **Q2:** Why do I get an **Illegal instruction (core dump)** error if I distribute PopIns2 jobs among a HPC cluster? **A:** PopIns2, by default, uses the `-march=native` compiler option to build the binary. This option opitimizes the code using the current maschine's processor architecture. If you distribute PopIns2 jobs among cluster nodes you have to make sure that all nodes support the same CPU instructions like the machine that the binary was built on.
