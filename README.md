# PopIns2
Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

### Requirements:

- GCC (vers. 5.4.1)
- [SeqAn](https://www.seqan.de/) (vers. 2.3.2)
- [Bifrost](https://github.com/pmelsted/bfgraph)

### Build:

```
git clone https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

### Build Documentation (optional):

You can auto-build a _html_ and _LaTeX_ code documentation with doxygen by
```
cd doc
doxygen Doxyfile
```

### Run:
( Work in progress)
( important reference) [KmerStream](https://github.com/pmelsted/KmerStream)

### Help:

```
./popins2 -h

Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs
================================================================

SYNOPSIS
    ./popins2 COMMAND [OPTIONS]

COMMAND
    single          Build a single sample compacted de Bruijn Graph.
    merge           Merge many samples into a colored compacted de Bruijn Graph.

VERSION
    0.3.0-9d07d9e, Date: on 2018-04-10 11:10:51

Try `./popins2 COMMAND --help' for more information on each command.
```

