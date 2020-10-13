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

Prior to the installation make sure your system meets all the requirements. For the default settings of PopIns2 a *Bifrost* installation with MAX_KMER_SIZE=64 is required. If the executables of the software dependencies (bwa, samtools, sickle) are not accessible systemwide, you have to write the full paths to the executables into a configfile (see Installation). Submodules come by default with the git clone, there is no need for a manual installation. For backward compatibility PopIns2 still offers to use the *Velvet assembler* (see [popins](https://github.com/bkehr/popins) for installation recommendation).

## Installation:

```
git clone --recursive https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

If you don't append the binaries of the software dependencies to your system `PATH`, set the paths to the binaries within the *popins2.config* prior to executing `make`. After the compilation with `make` you should see the binary *popins2* in the main folder.

## Usage:

PopIns2 is a tool consisting of several submodules. To display the help page of every module type _popins2 [command] --help_ as shown in the [help section](#help). 

#### The assemble command
```
popins2 assemble [OPTIONS] myMappedSample.bam
```
The assemble command filters, clips and assembles a set of mapped reads. At first, the unmapped reads from the input BAM file are identified.
The unmapped reads can be further filtered via an additional reference (parameter __-r__), for instance to remove contamination from the read library.
Next, if the option is specified, sequencing adapters will be removed from the remaining unmapped reads.
Finally, the unmapped reads are assembled using the minia multi-k pipeline (default) or the Velvet assembler. \
After completion the assemble module returns the assembly (called 'assembly_final.contigs.fa') and some meta files.

#### The merge command
```
popins2 merge [OPTIONS] {-s|-r} /path/to/input/
```
The merge command builds and writes a colored and compacted de Bruijn Graph (ccdbg) of all samples in a given source directory.
The source files in that directory can be of any (gzipped) plain sequence format: fa, fa.gz, fq, fq.gz, fasta, fasta.gz, fastq or fastq.gz.
Once the ccdbg is built, the merge module finds paths in the graph and returns them as contigs. \
It is common practise to create a directory of symbolic links to the assemblies from the _assemble_ module and use that directory as input.
In that case you probably like to use the __-r__ option since contigs out of an assembler are typically of high quality.
If you use raw reads for the merge module you probably like to use the __-s__ option to filter for low frequent k-mers. \
After completion the merge module returns the ccdbg and a file of merged contigs from all samples (called 'supercontigs').

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

- **Q:** Where do I install _SeqAn_ and _Bifrost_ for PopIns2? **A:** The C++ libraries SeqAn and Bifrost have to be found by your compiler, i.e. they should usually be located in your system's `local/include` and `local/lib` folders. The respective website and github page (see Requirements) provide more details.
