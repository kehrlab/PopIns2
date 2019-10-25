# PopIns2
Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

## Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 16.04 LTS, CentOS Linux 7.6 |
| C++14 capable compiler | g++ vers. 4.9.2, 5.5.0 |
| [SeqAn](https://www.seqan.de/) | vers. 2.3.2, 2.4.1 |
| [Bifrost](https://github.com/pmelsted/bfgraph) | vers. 0.2-699bdcf |

## Build:

```
git clone https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

## Usage:

PopIns2 is a tool consisting of many submodules. To display the help page of every module type _popins2 [command] --help_ as shown in the [help section](#help). 

#### The merge command
```
popins2 merge {-s|-r} /path/to/input/
```
The merge command builds and writes a colored and compacted de Bruijn Graph (ccdbg) of all samples in a given source directory. The source files in that directory can be of any (gzipped) plain sequence format: fa, fa.gz, fq, fq.gz, fasta, fasta.gz, fastq or fastq.gz. Once the ccdbg is built, the merge module finds paths in the graph and returns them as contigs. \
It is common practise to create a directory of symbolic links to the assemblies from the _assembly_ module and use that directory as input. In that case you probably like to use the __-r__ option since contigs out of an assembler are typically of high quality. If you use raw reads for the merge module you probably like to use the __-s__ option to filter for low frequent k-mers. \
After completion the merge module returns the ccdbg and a file of merged contigs from all samples (called 'supercontigs').

## Help:

```
./popins2 -h

Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs
================================================================

SYNOPSIS
    popins2 COMMAND [OPTIONS]

COMMAND
    merge           Merge many samples into a colored compacted de Bruijn Graph.

VERSION
    0.6.0-b7aee98, Date: on 2019-10-24 16:03:33

Try `popins2 COMMAND --help' for more information on each command.

```