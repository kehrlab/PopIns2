/**
 * @file    src/argument_parsing.h
 * @brief   Library for the command line parsing, user interaction and directory navigation.
 *
 */
#ifndef ARGUMENT_PARSING_H_
#define ARGUMENT_PARSING_H_


#include "util.h"
#include <bifrost/ColoredCDBG.hpp>
#include <seqan/arg_parse.h>
#include <vector>

using namespace std;



// =========================
// Option wrapper classes
// =========================

struct AssemblyOptions {
    CharString mappingFile;
    CharString matepairFile;
    CharString referenceFile;

    CharString prefix;
    CharString sampleID;

    unsigned kmerLength;
    CharString adapters;
    int humanSeqs;

    unsigned threads;
    CharString memory;

    bool use_velvet;
    bool use_spades;
    bool skip_assembly;
    float alignment_score_factor;

    AssemblyOptions () :
        matepairFile(""),
        referenceFile(""),
        prefix("."),
        sampleID(""),
        kmerLength(47),
        humanSeqs(maxValue<int>()),
        threads(1),
        memory("768M"),
        use_velvet(false),
        use_spades(false),
        skip_assembly(false),
        alignment_score_factor(0.67f)
    {}
};


struct MergeOptions {
    /************
    *  Bifrost  *
    ************/
    bool verbose;

    size_t nb_threads;

    vector<string> filename_seq_in;
    vector<string> filename_ref_in;

    std::string filename_graph_in;
    std::string filename_colors_in;

    int k;
    int g;

    bool clipTips;
    bool deleteIsolated;
    bool useMercyKmers;

    std::string prefixFilenameOut;
    std::string contigsFileName;

    /************
    *  PopIns2  *
    ************/
    int setcover_min_kmers;
    float min_entropy;
    bool write_setcover;
    bool write_lecc;

    MergeOptions () :       // the initializer list defines the program defaults
        verbose(false),
        nb_threads(1),

        filename_graph_in(""),
        filename_colors_in(""),

        k(63),
        g(49),
        clipTips(false),
        deleteIsolated(false),
        useMercyKmers(false),

        prefixFilenameOut("supercontigs"),
        contigsFileName("assembly_final.contigs.fa"),

        setcover_min_kmers(62),
        min_entropy(0.0f),
        write_setcover(false),
        write_lecc(false)
    {}
};


struct MultikOptions {
    int k_init;
    unsigned k_max;
    unsigned delta_k;

    vector<string> inputFiles;

    std::string samplePath;
    std::string tempPath;

    std::string prefixFilenameOut;

    MultikOptions () :
        k_init(27),
        k_max(127),
        delta_k(20),
        samplePath(""),
        tempPath("auxMultik"),
        prefixFilenameOut("ccdbg")
    {}
};


struct ContigMapOptions {
    CharString prefix;
    CharString sampleID;
    CharString contigFile;
    CharString referenceFile;

    bool bestAlignment;
    int maxInsertSize;
    bool deleteNonRefNew;

    unsigned threads;
    CharString memory;

    ContigMapOptions() :
        prefix("."), sampleID(""), contigFile("supercontigs.fa"), referenceFile("genome.fa"),
        bestAlignment(false), maxInsertSize(800), deleteNonRefNew(false), threads(1), memory("768M")
    {}
};


struct RefAlign_;
typedef Tag<RefAlign_> RefAlign;
struct SplitAlign_;
typedef Tag<SplitAlign_> SplitAlign;
struct SplitCombine_;
typedef Tag<SplitCombine_> SplitCombine;


template<typename TTag>
struct PlacingOptions {
    CharString prefix;
    CharString sampleID;
    CharString outFile;

    CharString locationsFile;
    CharString groupsFile;
    CharString supercontigFile;
    CharString referenceFile;

    double minLocScore;
    unsigned minAnchorReads;
    unsigned readLength;
    unsigned maxInsertSize;
    unsigned groupDist;

    PlacingOptions() :
        prefix("."), sampleID(""), outFile("insertions.vcf"), locationsFile("locations.txt"), groupsFile("groups.txt"),
        supercontigFile("supercontigs.fa"), referenceFile("genome.fa"),
        minLocScore(0.3), minAnchorReads(2), readLength(100), maxInsertSize(800), groupDist(100)
    {}
};


struct GenotypingOptions {
    CharString prefix;
    CharString sampleID;

    CharString referenceFile;
    CharString supercontigFile;
    CharString vcfFile;

    CharString genotypingModel;
    int regionWindowSize;
    bool addReadGroup;

    int maxInsertSize;
    int bpQclip;
    int minSeqLen;
    double minReadProb;
    int maxBARcount;

    int match;
    int mismatch;
    int gapOpen;
    int gapExtend;
    int minAlignScore;

    // hidden options
    bool verbose;
    bool callBoth;
    bool useReadCounts;
    bool fullOverlap;

    GenotypingOptions() :
        prefix("."), sampleID(""), referenceFile("genome.fa"), supercontigFile("supercontigs.fa"), vcfFile("insertions.vcf"),
      genotypingModel("RANDOM"), regionWindowSize(50), addReadGroup(false),
        maxInsertSize(500), bpQclip(0), minSeqLen(10), minReadProb(0.00001), maxBARcount(200),
        match(1), mismatch(-2), gapOpen(-4), gapExtend(-1), minAlignScore(55),
      verbose(false), callBoth(false), useReadCounts(false), fullOverlap(false)
    {}
};


// =========================
// Option transfer functions
// =========================

bool getOptionValues(AssemblyOptions & options, ArgumentParser const & parser){

    getArgumentValue(options.mappingFile, parser, 0);

    if (isSet(parser, "prefix"))
       getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "sample"))
       getOptionValue(options.sampleID, parser, "sample");
    if (isSet(parser, "matePair"))
        getOptionValue(options.matepairFile, parser, "matePair");
    if (isSet(parser, "adapters"))
        getOptionValue(options.adapters, parser, "adapters");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "filter"))
        getOptionValue(options.humanSeqs, parser, "filter");
    if (isSet(parser, "use-velvet"))
        getOptionValue(options.use_velvet, parser, "use-velvet");
    if (isSet(parser, "use-spades"))
        getOptionValue(options.use_spades, parser, "use-spades");
    if (isSet(parser, "skip-assembly"))
        getOptionValue(options.skip_assembly, parser, "skip-assembly");
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.kmerLength, parser, "kmerLength");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");
    if (isSet(parser, "alignment-score-factor"))
        getOptionValue(options.alignment_score_factor, parser, "alignment-score-factor");

    return true;
}


bool getOptionValues(MergeOptions &options, seqan::ArgumentParser &parser){

    if (isSet(parser, "verbose"))
        getOptionValue(options.verbose, parser, "verbose");
    if (isSet(parser, "threads"))
        getOptionValue(options.nb_threads, parser, "threads");
    if (isSet(parser, "contigs-filename"))
        getOptionValue(options.contigsFileName, parser, "contigs-filename");

    if (isSet(parser, "input-seq-files")){
        string indir;
        getOptionValue(indir, parser, "input-seq-files");
        listFiles(options.filename_seq_in, indir, options.contigsFileName);
    }
    if (isSet(parser, "input-ref-files")){
        string indir;
        getOptionValue(indir, parser, "input-ref-files");
        listFiles(options.filename_ref_in, indir, options.contigsFileName);
    }

    if (isSet(parser, "input-graph-file")){
        getOptionValue(options.filename_graph_in, parser, "input-graph-file");
    }
    if (isSet(parser, "input-colors-file")){
        getOptionValue(options.filename_colors_in, parser, "input-colors-file");
    }

    if (isSet(parser, "kmer-length"))
        getOptionValue(options.k, parser, "kmer-length");
    if (isSet(parser, "minimizer-length"))
        getOptionValue(options.g, parser, "minimizer-length");
    if (isSet(parser, "clip-tips"))
        getOptionValue(options.clipTips, parser, "clip-tips");
    if (isSet(parser, "del-isolated"))
        getOptionValue(options.deleteIsolated, parser, "del-isolated");
    if (isSet(parser, "mercy-kmers"))
        getOptionValue(options.useMercyKmers, parser, "mercy-kmers");
    if (isSet(parser, "outputfile-prefix"))
        getOptionValue(options.prefixFilenameOut, parser, "outputfile-prefix");

    if (isSet(parser, "setcover-min-kmers"))
        getOptionValue(options.setcover_min_kmers, parser, "setcover-min-kmers");
    if (isSet(parser, "min-entropy"))
        getOptionValue(options.min_entropy, parser, "min-entropy");
    if (isSet(parser, "write-setcover"))
        getOptionValue(options.write_setcover, parser, "write-setcover");
    if (isSet(parser, "write-lecc"))
        getOptionValue(options.write_lecc, parser, "write-lecc");

    return true;
}


/**
 *          This function transfers all relevant options of the merge module to
 *          a CCDBG_Build_opt instance.
 * @param   options is an options instance of the merge module
 * @param   graph_opt is an options instance of Bifrost
 */
void setupBifrostOptions(const MergeOptions &options, CCDBG_Build_opt &graph_opt){
    graph_opt.verbose            = options.verbose;
    graph_opt.nb_threads         = options.nb_threads;
    graph_opt.filename_seq_in    = options.filename_seq_in;
    graph_opt.filename_ref_in    = options.filename_ref_in;
    graph_opt.filename_graph_in  = options.filename_graph_in;
    graph_opt.filename_colors_in = options.filename_colors_in;
    graph_opt.k                  = options.k;
    graph_opt.g                  = options.g;
    graph_opt.clipTips           = options.clipTips;
    graph_opt.deleteIsolated     = options.deleteIsolated;
    graph_opt.useMercyKmers      = options.useMercyKmers;
    graph_opt.prefixFilenameOut  = options.prefixFilenameOut;
}


bool getOptionValues(MultikOptions &options, seqan::ArgumentParser &parser){
    if (isSet(parser, "sample-path")){
        // store sample path
        getOptionValue(options.samplePath, parser, "sample-path");
        // get fastx files from sample path
        getFastx(options.inputFiles, options.samplePath);
    }
    if (isSet(parser, "temp-path"))
        getOptionValue(options.tempPath, parser, "temp-path");
    if (isSet(parser, "outputfile-prefix"))
        getOptionValue(options.prefixFilenameOut, parser, "outputfile-prefix");
    if (isSet(parser, "k-init"))
        getOptionValue(options.k_init, parser, "k-init");
    if (isSet(parser, "k-max"))
        getOptionValue(options.k_max, parser, "k-max");
    if (isSet(parser, "delta-k"))
        getOptionValue(options.delta_k, parser, "delta-k");

    return true;
}


bool getOptionValues(ContigMapOptions &options, seqan::ArgumentParser &parser){
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.contigFile, parser, "contigs");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "best"))
        options.bestAlignment = true;
    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "noNonRefNew"))
        options.deleteNonRefNew = true;
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");

    return true;
}


bool getOptionValues(PlacingOptions<RefAlign> &options, seqan::ArgumentParser &parser){
    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "locations"))
        getOptionValue(options.locationsFile, parser, "locations");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "insertions"))
        getOptionValue(options.outFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "groups"))
        getOptionValue(options.groupsFile, parser, "groups");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "readLength"))
        getOptionValue(options.readLength, parser, "readLength");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minLocScore, parser, "minScore");
    if (isSet(parser, "minReads"))
        getOptionValue(options.minAnchorReads, parser, "minReads");
    if (isSet(parser, "groupDist"))
        getOptionValue(options.groupDist, parser, "groupDist");

    return true;
}


bool getOptionValues(PlacingOptions<SplitAlign> &options, seqan::ArgumentParser &parser){
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "readLength"))
        getOptionValue(options.readLength, parser, "readLength");

    return true;
}


bool getOptionValues(PlacingOptions<SplitCombine> &options, seqan::ArgumentParser &parser){
    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "insertions"))
        getOptionValue(options.outFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");

    return true;
}


bool getOptionValues(GenotypingOptions &options, seqan::ArgumentParser &parser){
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "insertions"))
        getOptionValue(options.vcfFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");

    if (isSet(parser, "model"))
        getOptionValue(options.genotypingModel, parser, "model");
    if(isSet(parser, "window"))
        getOptionValue( options.regionWindowSize, parser, "window");
    options.addReadGroup = isSet(parser, "addReadGroup");

    if (isSet(parser, "match"))
        getOptionValue(options.match, parser, "match");
    if (isSet(parser, "mismatch"))
        getOptionValue(options.mismatch, parser, "mismatch");
    if (isSet(parser, "gapOpen"))
        getOptionValue(options.gapOpen, parser, "gapOpen");
    if (isSet(parser, "gapExtend"))
        getOptionValue(options.gapExtend, parser, "gapExtend");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minAlignScore, parser, "minScore");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "minReadProb"))
        getOptionValue(options.minReadProb, parser, "minReadProb");
    if (isSet(parser, "maxReadCount"))
        getOptionValue(options.maxBARcount, parser, "maxReadCount");
    if (isSet(parser, "qual"))
        getOptionValue(options.bpQclip, parser, "qual");
    if (isSet(parser, "minSeqLen"))
        getOptionValue(options.minSeqLen, parser, "minSeqLen");

    options.verbose = isSet(parser, "verbose");
    options.callBoth = isSet(parser, "callBoth");
    options.fullOverlap = isSet(parser, "fullOverlap");
    options.useReadCounts = isSet(parser, "readCounts");

    return true;
}


// =========================
// Hide options functions
// =========================

void setHiddenOptions(ArgumentParser & parser, bool hide, AssemblyOptions &){
    hideOption(parser, "matePair", hide);
    hideOption(parser, "kmerLength", hide);
    hideOption(parser, "alignment-score-factor", hide);
}


void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, MergeOptions &){
    hideOption(parser, "minimizer-length",   hide);
    hideOption(parser, "mercy-kmers",        hide);

    hideOption(parser, "setcover-min-kmers", hide);
    hideOption(parser, "write-setcover",     hide);
    hideOption(parser, "write-lecc",         hide);
}


void setHiddenOptions(seqan::ArgumentParser & /*parser*/, bool /*hide*/, MultikOptions &){
    // TODO
}


void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, ContigMapOptions &){
   hideOption(parser, "b", hide);
   hideOption(parser, "e", hide);
   hideOption(parser, "d", hide);
}


void setHiddenOptions(seqan::ArgumentParser & parser, bool hide, PlacingOptions<RefAlign> &){
   hideOption(parser, "groupDist", hide);
}

void setHiddenOptions(seqan::ArgumentParser & /*parser*/, bool /*hide*/, PlacingOptions<SplitAlign> &){
	// Nothing to be done.
}

void setHiddenOptions(seqan::ArgumentParser & /*parser*/, bool /*hide*/, PlacingOptions<SplitCombine> &){
	// Nothing to be done.
}


void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, GenotypingOptions &){
   hideOption(parser, "addReadGroup", hide);

   hideOption(parser, "maxInsertSize", hide);
   hideOption(parser, "qual", hide);
   hideOption(parser, "minSeqLen", hide);
   hideOption(parser, "minReadProb", hide);
   hideOption(parser, "maxReadCount", hide);

   hideOption(parser, "match", hide);
   hideOption(parser, "mismatch", hide);
   hideOption(parser, "gapOpen", hide);
   hideOption(parser, "gapExtend", hide);
   hideOption(parser, "minScore", hide);
}


// ==========================================================================
// Functions setupParser()
// ==========================================================================

void setupParser(ArgumentParser & parser, AssemblyOptions & options){
    setShortDescription(parser, "Assembly of unmapped reads.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM_FILE\\fP");
    addDescription(parser, "Finds reads without high-quality alignment in the \\fIBAM FILE\\fP, quality filters them "
          "using SICKLE and assembles them into contigs using MINIA (default) or VELVET. If the option "
          "\'--reference \\fIFASTA FILE\\fP\' is set, the reads are first remapped to this reference using BWA-MEM and "
          "only reads that remain without high-quality alignment after remapping are quality-filtered and assembled.");

    // Require a bam file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("s", "sample", "An ID for the sample.", ArgParseArgument::STRING, "SAMPLE_ID"));
    addOption(parser, ArgParseOption("mp", "matePair", "(Currently only available for Velvet.)", ArgParseArgument::INPUT_FILE, "BAM FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("a", "adapters", "Enable adapter removal for Illumina reads. Default: \\fIno adapter removal\\fP.", ArgParseArgument::STRING, "STR"));
    addOption(parser, ArgParseOption("r", "reference", "Remap reads to this reference before assembly. Default: \\fIno remapping\\fP.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("f", "filter", "Treat reads aligned to all but the first INT reference sequences after remapping as high-quality aligned even if their alignment quality is low. "
          "Recommended for non-human reference sequences.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("vel", "use-velvet", "Use the velvet assembler. Default: Minia."));
    addOption(parser, ArgParseOption("d", "use-spades", "Use the SPAdes assembler. Default: Minia."));
    addOption(parser, ArgParseOption("n", "skip-assembly", "Skip assembly per sample."));
    addOption(parser, ArgParseOption("k", "kmerLength", "The k-mer size if the velvet assembler is used.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("c", "alignment-score-factor", "A record is considered low quality if the alignment score (AS) is below FLOAT*read length", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    addSection(parser, "Compute resource options");
    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for BWA and samtools sort.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory per thread for samtools sort; suffix K/M/G recognized.", ArgParseArgument::STRING, "STR"));

    // Set valid and default values.
    setValidValues(parser, "adapters", "HiSeq HiSeqX");
    setValidValues(parser, "reference", "fa fna fasta gz");
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "sample", "retrieval from BAM file header");
    setDefaultValue(parser, "kmerLength", options.kmerLength);
    setDefaultValue(parser, "threads", options.threads);
    setDefaultValue(parser, "memory", options.memory);
    setDefaultValue(parser, "alignment-score-factor", options.alignment_score_factor);

    setMinValue(parser, "threads", "1");
    setMinValue(parser, "alignment-score-factor", "0.0");
    setMaxValue(parser, "alignment-score-factor", "1.0");

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


/**
 *          Function to handle the input parsing for the popins2 merge module
 * @param   parser is a seqan argument parser instance
 * @param   options is a struct to store the input arguments for the merge module
 */
void setupParser(seqan::ArgumentParser &parser, MergeOptions &options){
    // Setup meta-information
    seqan::setShortDescription(parser, "Build or read a colored and compacted de Bruijn Graph (CCDBG) and generate supercontigs.");
    seqan::setVersion(parser, VERSION);
    seqan::setDate(parser, DATE);
    seqan::addUsageLine(parser, "\\--input-{seq|ref}-files DIR or --input-graph-file GFA --input-colors-file BFG_COLORS [OPTIONS] \\fP ");

    // Setup options
    seqan::addSection(parser, "I/O options");
    seqan::addOption(parser, seqan::ArgParseOption("s", "input-seq-files",   "Path to the sample directories.", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("r", "input-ref-files",   "Path to the sample directories. (no abundance filter)", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("y", "input-graph-file",  "Source file with dBG", seqan::ArgParseArgument::STRING, "GFA"));
    seqan::addOption(parser, seqan::ArgParseOption("z", "input-colors-file", "Source file with dBG colors", seqan::ArgParseArgument::STRING, "BFG_COLORS"));
    seqan::addOption(parser, seqan::ArgParseOption("p", "outputfile-prefix", "Specify a prefix for the output files.", seqan::ArgParseArgument::STRING, "STRING"));
    seqan::addOption(parser, seqan::ArgParseOption("f", "contigs-filename",  "Specify a filename of contigs to search for in the sample directories.", seqan::ArgParseArgument::STRING, "STRING"));
    seqan::addOption(parser, seqan::ArgParseOption("c", "write-setcover",    "Write a CSV file with unitig IDs of the setcover"));
    seqan::addOption(parser, seqan::ArgParseOption("l", "write-lecc",        "Write a CSV file with unitig IDs of the LECCs"));

    seqan::addSection(parser, "Algorithm options");
    seqan::addOption(parser, seqan::ArgParseOption("k", "kmer-length",        "Kmer length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("g", "minimizer-length",   "Minimizer length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("v", "verbose",            "Print more output"));
    seqan::addOption(parser, seqan::ArgParseOption("i", "clip-tips",          "Clip tips shorter than k kmers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("d", "del-isolated",       "Delete isolated contigs shorter than k kmers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("x", "mercy-kmers",        "Keep low coverage k-mers (cov=1) connecting tips of the graph"));
    seqan::addOption(parser, seqan::ArgParseOption("m", "setcover-min-kmers", "Minimum amount of unseen kmers to include a path into the set cover", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("e", "min-entropy",        "Minimum entropy for a unitig to not get flagged as low entropy", seqan::ArgParseArgument::DOUBLE, "FLOAT"));

    seqan::addSection(parser, "Compute resource options");
    seqan::addOption(parser, seqan::ArgParseOption("t", "threads", "Amount of threads for parallel processing", seqan::ArgParseArgument::INTEGER, "INT"));

    // Setup option constraints
    seqan::setDefaultValue(parser, "p",   options.prefixFilenameOut);
    seqan::setDefaultValue(parser, "f",   options.contigsFileName);
    seqan::setDefaultValue(parser, "k",   options.k);
    seqan::setDefaultValue(parser, "g",   options.g);
    seqan::setDefaultValue(parser, "m",   options.setcover_min_kmers);
    seqan::setDefaultValue(parser, "e",   options.min_entropy);
    seqan::setDefaultValue(parser, "t",   options.nb_threads);

    seqan::setMinValue(parser, "k", "1");
    seqan::setMaxValue(parser, "k", std::to_string(MAX_KMER_SIZE-1));
    seqan::setMinValue(parser, "g", "1");
    seqan::setMaxValue(parser, "g", std::to_string(MAX_KMER_SIZE-3));
    seqan::setMinValue(parser, "m", "1");
    seqan::setMinValue(parser, "e", "0.0");
    seqan::setMaxValue(parser, "e", "1.0");
    seqan::setMinValue(parser, "t", "1");

    // Setup hidden options
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, MultikOptions &options){
    // Setup meta-information
    seqan::setShortDescription(parser, "Multi-k framework for a colored and compacted de Bruijn Graph (CCDBG)");
    seqan::setVersion(parser, VERSION);
    seqan::setDate(parser, DATE);
    seqan::addUsageLine(parser, "\\--sample-path STRING [OPTIONS] \\fP ");

    // Setup options
    seqan::addSection(parser, "I/O options");
    seqan::addOption(parser, seqan::ArgParseOption("s", "sample-path",   "Source directory with FASTA/Q files", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("a", "temp-path",     "Auxiliary directory for temporary files.", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("p", "outputfile-prefix", "Specify a prefix for the output files", seqan::ArgParseArgument::STRING, "STRING"));

    seqan::addSection(parser, "Algorithm options");
    seqan::addOption(parser, seqan::ArgParseOption("k", "k-init",    "Initial kmer length to start the multi-k iteration", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("m", "k-max",     "Maximal kmer length to build a dBG with", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("d", "delta-k",   "Step size to increase k", seqan::ArgParseArgument::INTEGER, "INT"));

    // Setup option constraints
    seqan::setDefaultValue(parser, "p",   options.prefixFilenameOut);
    seqan::setDefaultValue(parser, "a",   options.tempPath);
    seqan::setDefaultValue(parser, "k",   options.k_init);
    seqan::setDefaultValue(parser, "m",   options.k_max);
    seqan::setDefaultValue(parser, "d",   options.delta_k);

    // Setup hidden options
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, ContigMapOptions &options){
    setShortDescription(parser, "Alignment of unmapped reads to assembled contigs.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");
    addDescription(parser, "Aligns the reads with low-quality alignments of a sample to the set of supercontigs using "
            "BWA-MEM. Merges the BWA output file with the sample's non_ref.bam file into a non_ref_new.bam file where "
    		"information about read mates is set.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of (super-)contigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("b", "best", "Do not use BWA-mem's -a option to output all alignments of a read."));
    addOption(parser, ArgParseOption("e", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("d", "noNonRefNew", "Delete the non_ref_new.bam file after writing locations."));

    addSection(parser, "Compute resource options");
    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for BWA and samtools sort.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory per thread for samtools sort; suffix K/M/G recognized.", ArgParseArgument::STRING, "STR"));

    // Set valid values.
    setMinValue(parser, "threads", "1");
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "contigs", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "contigs", options.contigFile);
    setDefaultValue(parser, "reference", options.referenceFile);
    setDefaultValue(parser, "best", "false");
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);
    setDefaultValue(parser, "noNonRefNew", "false");
    setDefaultValue(parser, "threads", options.threads);
    setDefaultValue(parser, "memory", options.memory);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, PlacingOptions<RefAlign> &options){
    setShortDescription(parser, "Contig placing by alignment of contig ends to reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP]");

    addDescription(parser, "This is step 1/3 of contig placing. The contig locations in the sample directories are "
    		"merged into one file of locations. Next, prefixes/suffixes of contigs are aligned to the merged locations "
    		"on the reference genome and VCF records are written if the alignment is successful. Locations of contigs "
    		"that do not align to the reference genome are written to additional output files \\fIlocations_unplaced.txt\\fP "
    		"in the sample directories. Further, groups of contigs that can be placed at the same position and whose "
    		"prefixes/suffixes align to each other are written to another output file; only a single VCF record is "
    		"written per group.");

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("l", "locations", "Name of merged locations file.", ArgParseArgument::OUTPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF output file.", ArgParseArgument::OUTPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("g", "groups", "Name of file containing groups of contigs that can be placed at the same position and whose prefixes/suffixes align.", ArgParseArgument::OUTPUT_FILE, "FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("", "minScore", "Minimal anchoring score for a location.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("", "minReads", "Minimal number of anchoring read pairs for a location.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "groupDist", "Minimal distance between groups of locations.", ArgParseArgument::INTEGER, "INT"));

    // Set valid values.
    setMinValue(parser, "minScore", "0");
    setMaxValue(parser, "minScore", "1");
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "locations", options.locationsFile);
    setDefaultValue(parser, "insertions", options.outFile);
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);
    setDefaultValue(parser, "groups", options.groupsFile);

    setDefaultValue(parser, "minScore", options.minLocScore);
    setDefaultValue(parser, "minReads", options.minAnchorReads);
    setDefaultValue(parser, "groupDist", options.groupDist);
    setDefaultValue(parser, "readLength", options.readLength);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, PlacingOptions<SplitAlign> &options){
    setShortDescription(parser, "Contig placing by split-read alignment.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");

    addDescription(parser, "This is step 2/3 of contig placing. All locations in a sample's "
    		"\\fIlocations_unplaced.txt\\fP are split-read aligned and the results are written to a file "
    		"\\fIlocations_placed.txt\\fP in the sample directory.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));

    // Set valid values.
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    setDefaultValue(parser, "readLength", options.readLength);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, PlacingOptions<SplitCombine> &options){
    setShortDescription(parser, "Combining breakpoint positions from split-read alignment.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP]");

    addDescription(parser, "This is step 3/3 of contig placing. The results from split-read alignment (the "
    		"\\fIlocations_placed.txt\\fP files) of all samples are combined and appended to the VCF output file.");

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF output file.", ArgParseArgument::OUTPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    // Set valid values.
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "insertions", options.outFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


void setupParser(seqan::ArgumentParser &parser, GenotypingOptions &options){
    setShortDescription(parser, "Genotyping a sample for insertions.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");
    addDescription(parser, "Computes genotype likelihoods for a sample for all insertions given in the input VCF file "
    		"by aligning all reads, which are mapped to the reference genome around the insertion breakpoint or to the "
    		"contig, to the reference and to the alternative insertion sequence. VCF records with the genotype "
    		"likelihoods in GT:PL format for the individual are written to a file \\fIinsertions.vcf\\fP in the sample "
            "directory.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF input file.", ArgParseArgument::INPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("m", "model", "Model used for genotyping.", ArgParseArgument::STRING, "GENOTYPING_MODEL"));
    addOption(parser, ArgParseOption("w", "window", "Region window size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("rg", "addReadGroup", "Add read group."));

    addSection(parser, "Read(-pair) options");
    addOption(parser, ArgParseOption("", "maxInsertSize", "Maximum read pair insert size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "qual", "Quality score threshold for read trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minSeqLen", "Minimum read length after trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minReadProb", "Minimum read probability.", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("", "maxReadCount", "Maximum number of reads to consider in region window.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Alignment options");
    addOption(parser, ArgParseOption("", "match", "Cost of match.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "mismatch", "Cost of mismatch.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "gapOpen", "Cost of gap open.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "gapExtend", "Cost of gap extend.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minScore", "Minimum alignment score.", ArgParseArgument::INTEGER, "INT"));

    // Misc hidden options.
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
   hideOption(parser, "verbose", true);
    addOption(parser, ArgParseOption("", "callBoth", "Call both models."));
   hideOption(parser, "callBoth", true);
    addOption(parser, ArgParseOption("", "readCounts", "Use read counts."));
   hideOption(parser, "readCounts", true);
    addOption(parser, ArgParseOption("", "fullOverlap", "Full overlap of read."));
   hideOption(parser, "fullOverlap", true);

    // Set valid values.
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");
    setValidValues(parser, "model", "DUP RANDOM");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "insertions", options.vcfFile);
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    setDefaultValue(parser, "model", options.genotypingModel);
    setDefaultValue(parser, "window", options.regionWindowSize);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);
    setDefaultValue(parser, "qual", options.bpQclip);
    setDefaultValue(parser, "minReadProb", options.minReadProb);
    setDefaultValue(parser, "maxReadCount", options.maxBARcount);

    setDefaultValue(parser, "match", options.match);
    setDefaultValue(parser, "mismatch", options.mismatch);
    setDefaultValue(parser, "gapOpen", options.gapOpen);
    setDefaultValue(parser, "gapExtend", options.gapExtend);
    setDefaultValue(parser, "minScore", options.minAlignScore);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


// ==========================================================================
// Function checkInput()
// ==========================================================================

ArgumentParser::ParseResult checkInput(AssemblyOptions & options){

    ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

    if (options.use_spades && options.use_velvet)
    {
        std::cerr << "ERROR: Multiple assemblers specified! Choose only one of them per program call." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.prefix != "." && !exists(options.prefix))
    {
        std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (!exists(options.mappingFile))
    {
        std::cerr << "ERROR: Input BAM file \'" << options.mappingFile << "\' does not exist." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    CharString baiFile = options.mappingFile;
    baiFile += ".bai";
    if (!exists(baiFile))
    {
        std::cerr << "ERROR: BAM index file \'" << baiFile << "\' does not exist." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.matepairFile != "" && !exists(options.matepairFile))
    {
        std::cerr << "ERROR: Input BAM file \'" << options.matepairFile << "\' does not exist." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    return res;
}


ArgumentParser::ParseResult checkInput(MergeOptions & options){

    ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

    if (options.filename_ref_in.empty() && options.filename_seq_in.empty() && strcmp(options.filename_graph_in.c_str(), "")==0 ){
        cerr << "[popins2 merge][parser] ERROR: No input files or graph found." << endl;
        cerr << "[popins2 merge][parser] ERROR: At least one input (-r/-s) must be given OR a graph (-y) AND colors (-z)." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if ((!options.filename_ref_in.empty() || !options.filename_seq_in.empty()) && strcmp(options.filename_graph_in.c_str(), "")!=0 ){
        cerr << "[popins2 merge][parser] ERROR: The merge module takes only sequence files or a graph, not both!" << endl;
        cerr << "[popins2 merge][parser] ERROR: At least one input (-r/-s) must be given XOR a graph (-y) AND colors (-z)." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (strcmp(options.filename_graph_in.c_str(), "")!=0 && strcmp(options.filename_colors_in.c_str(), "")==0  ||
        strcmp(options.filename_graph_in.c_str(), "")==0 && strcmp(options.filename_colors_in.c_str(), "")!=0){
        cerr << "[popins2 merge][parser] ERROR: One of the colored de Bruijn Graph files is missing (-y/-z)." << endl;
        cerr << "[popins2 merge][parser] ERROR: If a graph should be read, please provide a graph (-y) AND a colors file (-z)." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.g > options.k - 2){
        cerr << "[popins2 merge][parser] ERROR: Minimizer length (-g) should not be larger than k-2." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    return res;
}


ArgumentParser::ParseResult checkInput(MultikOptions & options){

    ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

    if (options.samplePath == "") {
        cerr << "[popins2 multik][parser] ERROR: No input path specified." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.k_init < 1) {
        cerr << "[popins2 multik][parser] ERROR: Parameter k_init must be a positive integer." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if ((unsigned)options.k_init > options.k_max) {
        cerr << "[popins2 multik][parser] ERROR: Parameter k_init must be smaller than k_max." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.delta_k > options.k_max) {
        cerr << "[popins2 multik][parser] ERROR: Parameter delta_k must be smaller than k_max." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    if (options.delta_k > (options.k_max - (unsigned)options.k_init)) {
        cerr << "[popins2 multik][parser] ERROR: This step size delta_k will have no effect." << endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    return res;
}


ArgumentParser::ParseResult checkInput(ContigMapOptions &options){

	ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

	if (options.prefix != "." && !exists(options.prefix))
	{
		std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.contigFile))
	{
		std::cerr << "ERROR: Contig file \'" << options.contigFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.referenceFile))
	{
		std::cerr << "ERROR: Reference genome file \'" << options.referenceFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	return res;
}


ArgumentParser::ParseResult checkInput(PlacingOptions<RefAlign> &options){

	ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

	if (options.prefix != "." && !exists(options.prefix))
	{
		std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.supercontigFile))
	{
		std::cerr << "ERROR: Contig file \'" << options.supercontigFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.referenceFile))
	{
		std::cerr << "ERROR: Reference genome file \'" << options.referenceFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	return res;
}


ArgumentParser::ParseResult checkInput(PlacingOptions<SplitAlign> &options){

	ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

	if (options.prefix != "." && !exists(options.prefix))
	{
		std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.supercontigFile))
	{
		std::cerr << "ERROR: Contig file \'" << options.supercontigFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.referenceFile))
	{
		std::cerr << "ERROR: Reference genome file \'" << options.referenceFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	return res;
}


ArgumentParser::ParseResult checkInput(PlacingOptions<SplitCombine> &options){

	ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

	if (options.prefix != "." && !exists(options.prefix))
	{
		std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.referenceFile))
	{
		std::cerr << "ERROR: Reference genome file \'" << options.referenceFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	return res;
}


ArgumentParser::ParseResult checkInput(GenotypingOptions &options){

	ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

	if (options.prefix != "." && !exists(options.prefix))
	{
		std::cerr << "ERROR: Path to sample directories \'" << options.prefix << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.vcfFile))
	{
		std::cerr << "ERROR: Input VCF file \'" << options.vcfFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.supercontigFile))
	{
		std::cerr << "ERROR: Contig file \'" << options.supercontigFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	if (!exists(options.referenceFile))
	{
		std::cerr << "ERROR: Reference genome file \'" << options.referenceFile << "\' does not exist." << std::endl;
		res = ArgumentParser::PARSE_ERROR;
	}

	return res;
}


// =========================
// Print functions
// =========================

void printHelp(char const * name){
    std::cerr << "Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs" << std::endl;
    std::cerr << "================================================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " COMMAND\033[0m [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mCOMMAND\033[0m" << std::endl;
    std::cerr << "    \033[1massemble\033[0m            Filter, clip and assemble unmapped reads from a sample." << std::endl;
    std::cerr << "    \033[1mmerge\033[0m               Generate supercontigs from a colored compacted de Bruijn Graph." << std::endl;
    std::cerr << "    \033[1mmultik\033[0m              Multi-k framework for a colored compacted de Bruijn Graph." << std::endl;
    std::cerr << "    \033[1mcontigmap\033[0m           Map unmapped reads to (super-)contigs." << std::endl;
    std::cerr << "    \033[1mplace-refalign\033[0m      Find position of (super-)contigs by aligning contig ends to the reference genome." << std::endl;
    std::cerr << "    \033[1mplace-splitalign\033[0m    Find position of (super-)contigs by split-read alignment (per sample)." << std::endl;
    std::cerr << "    \033[1mplace-finish\033[0m        Combine position found by split-read alignment from all samples." << std::endl;
    std::cerr << "    \033[1mgenotype\033[0m            Determine genotypes of all insertions in a sample." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << VERSION << ", Date: " << DATE << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}


void printMergeOptions(const MergeOptions &options){
    cout << "=========================================================" << endl;
    cout << "popins2 version    : " << VERSION                          << endl;
    cout << "PARAMETER ======== : VALUE ==============================" << endl;
    cout << "verbose            : " << options.verbose                  << endl;
    cout << "threads            : " << options.nb_threads               << endl;
    cout << "#filename_seq_in   : " << options.filename_seq_in.size()   << endl;
    cout << "#filename_ref_in   : " << options.filename_ref_in.size()   << endl;
    cout << "filename_graph_in  : " << options.filename_graph_in        << endl;
    cout << "filename_colors_in : " << options.filename_colors_in       << endl;
    cout << "k                  : " << options.k                        << endl;
    cout << "g                  : " << options.g                        << endl;
    cout << "clip-tips          : " << options.clipTips                 << endl;
    cout << "delete-isolated    : " << options.deleteIsolated           << endl;
    cout << "mercy-kmers        : " << options.useMercyKmers            << endl;
    cout << "outputfile-prefix  : " << options.prefixFilenameOut        << endl;
    cout << "setcover-min-kmers : " << options.setcover_min_kmers       << endl;
    cout << "min-entropy        : " << options.min_entropy              << endl;
    cout << "write-setcover     : " << options.write_setcover           << endl;
    cout << "write-lecc         : " << options.write_lecc               << endl;
    cout << "=========================================================" << endl;
}


void printMultikOptions(const MultikOptions &options){
    cout << "=========================================================" << endl;
    cout << "popins2 version    : " << VERSION                          << endl;
    cout << "PARAMETER ======== : VALUE ==============================" << endl;
    cout << "sample-path        : " << options.samplePath               << endl;
    cout << "#input files       : " << options.inputFiles.size()        << endl;
    cout << "temp-path          : " << options.tempPath                 << endl;
    cout << "outputfile-prefix  : " << options.prefixFilenameOut        << endl;
    cout << "k-init             : " << options.k_init                   << endl;
    cout << "k-max              : " << options.k_max                    << endl;
    cout << "delta-k            : " << options.delta_k                  << endl;
    cout << "=========================================================" << endl;
}


// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

/**
 *          Setup module based option forwarding
 * @param   options is an individual options struct instance, the type depends on the module called
 * @param   argc is the CLI argument count
 * @param   argv is the CLI argument list
 * @return  seqan::ArgumentParser::ParseResult type
 */
template<typename TOptions>
seqan::ArgumentParser::ParseResult parseCommandLine(TOptions &options, int argc, char const ** argv){

    // Concatenate program name from name and command.
    seqan::CharString prog_name = argv[0];
    prog_name += " ";
    prog_name += argv[1];

    ++argv;
    --argc;

    // Setup the parser.
    seqan::ArgumentParser parser(toCString(prog_name));
    seqan::addOption(parser, seqan::ArgParseOption("H", "fullHelp", "Display the help message with the full list of options."));
    setupParser(parser, options);

    // Parse the command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if full help option is set.
    if (isSet(parser, "fullHelp")){
        setHiddenOptions(parser, false, options);
        printHelp(parser);
        return seqan::ArgumentParser::PARSE_HELP;
    }

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Collect the option values.
    if (getOptionValues(options, parser) == false)
        return seqan::ArgumentParser::PARSE_HELP;

    // Check if input files exist.
    res = checkInput(options);

    return res;
}






#endif /*ARGUMENT_PARSING_H_*/
