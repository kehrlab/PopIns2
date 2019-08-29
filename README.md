# PopIns2
Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

### Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 16.04 LTS, CentOS Linux 7.6 |
| C++14 capable compiler | g++ vers. 4.9.2, 5.5.0 |
| [SeqAn](https://www.seqan.de/) | vers. 2.3.2 |
| [Bifrost](https://github.com/pmelsted/bfgraph) | vers. 0.2-699bdcf |

### Build:

```
git clone https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

### Run:

To build and write a colored compacted de Bruijn Graph you have to specify a source repository of (possibly gzipped) FAST[A/Q] files and an output prefix, e.g.:
```
./popins2 merge -s /path/to/input/ -o outPrefix -div
```

### Help:

```
./popins2 -h

Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs
================================================================

SYNOPSIS
    ./popins2 COMMAND [OPTIONS]

COMMAND
    merge           Merge many samples into a colored compacted de Bruijn Graph.

VERSION
    0.3.0-83b82e1, Date: on 2018-06-01 13:39:22

Try `./popins2 COMMAND --help' for more information on each command.
```

