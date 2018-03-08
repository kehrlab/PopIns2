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

You can auto-build a html and latex code documentation with doxygen by
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
    popins2 --indir DIR --prefix STRING --unique-kmers INT --non-unique-kmers INT [OPTIONS] 

DESCRIPTION
    Population-scale detection of novel-sequence insertions using de Bruijn Graphs

OPTIONS
    -h, --help
          Display the help message.
    --version
          Display version information.
    -i, --indir STRING
          Source directory with FASTQ file(s).
    -o, --prefix STRING
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
    -v, --verbose
          Print more output
    -c, --clip-tips
          Remove ends of the dBG with length <k
    -r, --remove-isolated
          Remove single contigs with length <k
    -f, --bloom-filter-file STRING
          Filename for a binary file storing the bloom filter data structure.

VERSION
    Last update: 
    PopIns2 version: 0.1.0-f25ec89
    SeqAn version: 2.3.2
