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

### Code Documentation (optional):

You can auto-build a _html_ and _LaTeX_ code documentation with doxygen by
```
cd doc
doxygen Doxyfile
```

### Run:

To build and write a compacted de Bruijn Graph all that is mandatory is a source repository of FASTX files and an output prefix:
```
./popins2 single -f /path/to/input/ -o myOutfile
```
However, the typical <code>popins2 single</code> call to prepare for <code>popins2 merge</code> is the 'diva' call, adding several graph option flags:
```
./popins2 single -f /path/to/input/ -o myOutfile -diva --threads 4 > single.log
```
The <code>-a</code> flag is causing the single module to output a FASTA instead of a GFA file. FASTAs of multiple single calls can then be merged into a colored compacted de Bruijn Graph (cdBG):
```
./popins2 merge -f /path/to/fastas/ -o myMergeOutfile --threads 4 --verbose > merge.log
```
The resulting cdBG is stored in GFA format and can be visualized by tools like [Bandage](https://github.com/rrwick/Bandage).

Advanced options tip:

For the _single_ module it is recommended to use [KmerStream](https://github.com/pmelsted/KmerStream) in advance if many samples will be processed. Popins2 internally calls KmerStream per sample if the (hidden) parameters <code>-n</code> and <code>-N</code> are not set. To speed up <code>popins2 single</code> for many samples, it is beneficial to call KmerStream on the largest sample and pass <code>-n=F0</code> and <code>-N=F0-f1</code> to all _single_ module runs.

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

