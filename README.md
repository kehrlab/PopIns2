# PopIns2
Population-scale detection of novel-sequence insertions using de Bruijn Graphs

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
./popins2 --help


popins2
=======

SYNOPSIS
    popins2 --file-dir DIR --output-file STRING [OPTIONS] 

DESCRIPTION
    Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs

OPTIONS
    -h, --help
          Display the help message.
    --version
          Display version information.
    -f, --file-dir STRING
          Source directory with input FASTX file(s).
    -F, --tmp-file-dir STRING
          Source directory with intermediate FASTA file(s) from the popins2 single module.
    -o, --output-file STRING
          Prefix for the GFA file
    -n, --unique-kmers INTEGER
          Amount of unique kmers. In range [1..inf].
    -N, --non-unique-kmers INTEGER
          Amount of non-unique kmers. In range [1..inf].
    -k, --kmer-length INTEGER
          K-mer length for the dBG construction. In range [1..63]. Default: 31.
    -g, --minimizer-length INTEGER
          Minimizer-length for the dBG construction. In range [1..62]. Default: 23.
    -t, --threads INTEGER
          Amount of threads for parallel processing. In range [1..inf]. Default: 1.
    -i, --clip-tips
          Remove branching ends (tips) of the dBG shorter than k k-mers in length
    -d, --del-isolated
          Remove single contigs shorter than k k-mers in length
    -v, --verbose
          Print more output

VERSION
    Last update: 
    popins2 version: 0.3.0-b9bff3c
    SeqAn version: 2.3.2

