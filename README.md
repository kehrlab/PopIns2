# PopIns2
Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

### Requirements:

- GCC (vers. 5.5.0)
- [SeqAn](https://www.seqan.de/) (vers. 2.3.2)
- [Bifrost](https://github.com/pmelsted/bfgraph) (vers. 0.2-d479e63)

### Build:

```
git clone https://github.com/kehrlab/PopIns2.git
cd PopIns2
mkdir build
make
```

### Code Documentation (optional):

You can auto-build a _html_ and _LaTeX_ code documentation with doxygen by
```
cd doc
doxygen Doxyfile
```

### Run:

To build and write a colored compacted de Bruijn Graph you have to specify a source repository of (possibly gzipped) FAST[A/Q] files and an output prefix:
```
./popins2 merge -s /path/to/input/ -o outPrefix
```

(OLD, NEED FIX, TODO)

However, the typical <code>popins2 single</code> call to prepare for <code>popins2 merge</code> is the 'diva' call, adding several graph option flags:
```
./popins2 single -f /path/to/input/ -o myOutfile -diva --threads 4 > single.log
```

(TO DO)

Advanced options tip:

[KmerStream](https://github.com/pmelsted/KmerStream)
To speed up <code>popins2 single</code> for many samples, it is beneficial to call KmerStream on the largest sample and pass <code>-n=F0</code> and <code>-N=F0-f1</code> to all _single_ module runs.

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

