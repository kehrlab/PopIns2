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

    AssemblyOptions () :
        matepairFile(""),
        referenceFile(""),
        prefix("."),
        sampleID(""),
        kmerLength(47),
        humanSeqs(maxValue<int>()),
        threads(1),
        memory("768M"),
        use_velvet(false)
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

    string prefixFilenameOut;

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
        g(23),
        clipTips(false),
        deleteIsolated(false),
        useMercyKmers(false),
        prefixFilenameOut("ccdbg"),

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

    std::string samplePath;
    std::string tempPath;

    MultikOptions () :
        k_init(27),
        k_max(127),
        delta_k(20),
        samplePath(""),
        tempPath("auxMultik")
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
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.kmerLength, parser, "kmerLength");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");

    return true;
}


bool getOptionValues(MergeOptions &options, seqan::ArgumentParser &parser){

    if (isSet(parser, "verbose"))
        getOptionValue(options.verbose, parser, "verbose");
    if (isSet(parser, "threads"))
        getOptionValue(options.nb_threads, parser, "threads");
    if (isSet(parser, "input-seq-files")){
        string indir;
        getOptionValue(indir, parser, "input-seq-files");
        getFastx(options.filename_seq_in, indir, options.verbose);
    }
    if (isSet(parser, "input-ref-files")){
        string indir;
        getOptionValue(indir, parser, "input-ref-files");
        getFastx(options.filename_ref_in, indir, options.verbose);
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

    if (isSet(parser, "sample-path"))
        getOptionValue(options.samplePath, parser, "sample-path");
    if (isSet(parser, "temp-path"))
        getOptionValue(options.tempPath, parser, "temp-path");
    if (isSet(parser, "k-init"))
        getOptionValue(options.k_init, parser, "k-init");
    if (isSet(parser, "k-max"))
        getOptionValue(options.k_max, parser, "k-max");
    if (isSet(parser, "delta-k"))
        getOptionValue(options.delta_k, parser, "delta-k");


    return true;
}


// =========================
// Hide options functions
// =========================

void setHiddenOptions(ArgumentParser & parser, bool hide, AssemblyOptions &){
    hideOption(parser, "matePair", hide);
    hideOption(parser, "kmerLength", hide);
}


void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, MergeOptions &){
    hideOption(parser, "minimizer-length",   hide);
    hideOption(parser, "mercy-kmers",        hide);

    hideOption(parser, "setcover-min-kmers", hide);
    hideOption(parser, "write-setcover",     hide);
    hideOption(parser, "write-lecc",     hide);
}

void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, MultikOptions &){
    // TODO
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
    addOption(parser, ArgParseOption("vel", "use-velvet", "Use the velvet assembler."));
    addOption(parser, ArgParseOption("k", "kmerLength", "The k-mer size if the velvet assembler is used.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Compute resource options");
    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for BWA and samtools sort.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory per thread for samtools sort; suffix K/M/G recognized.", ArgParseArgument::STRING, "STR"));

    // Set valid and default values.
    setValidValues(parser, "adapters", "HiSeq HiSeqX");
    setValidValues(parser, "reference", "fa fna fasta gz");
    setMinValue(parser, "threads", "1");

    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "sample", "retrieval from BAM file header");
    setDefaultValue(parser, "kmerLength", options.kmerLength);
    setDefaultValue(parser, "threads", options.threads);
    setDefaultValue(parser, "memory", options.memory);

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
    seqan::addOption(parser, seqan::ArgParseOption("s", "input-seq-files",   "Source directory with FASTA/Q files", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("r", "input-ref-files",   "Source directory with reference FASTA/Q files (no abundance filter)", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("y", "input-graph-file",  "Source file with dBG", seqan::ArgParseArgument::STRING, "GFA"));
    seqan::addOption(parser, seqan::ArgParseOption("z", "input-colors-file", "Source file with dBG colors", seqan::ArgParseArgument::STRING, "BFG_COLORS"));
    seqan::addOption(parser, seqan::ArgParseOption("p", "outputfile-prefix", "Specify a prefix for the output files", seqan::ArgParseArgument::STRING, "STRING"));
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
    seqan::setDefaultValue(parser, "k",   options.k);
    seqan::setDefaultValue(parser, "g",   options.g);
    seqan::setDefaultValue(parser, "m",   options.setcover_min_kmers);
    seqan::setDefaultValue(parser, "e",   options.min_entropy);
    seqan::setDefaultValue(parser, "t",   options.nb_threads);

    seqan::setMinValue(parser, "k", "1");
    seqan::setMaxValue(parser, "k", "63");
    seqan::setMinValue(parser, "g", "1");
    seqan::setMaxValue(parser, "g", "62");
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

    seqan::addSection(parser, "Algorithm options");
    seqan::addOption(parser, seqan::ArgParseOption("k", "k-init",    "Initial kmer length to start the multi-k iteration", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("m", "k-max",     "Maximal kmer length to build a dBG with", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("d", "delta-k",   "Step size to increase k", seqan::ArgParseArgument::INTEGER, "INT"));

    // Setup option constraints
    seqan::setDefaultValue(parser, "a",   options.tempPath);
    seqan::setDefaultValue(parser, "k",   options.k_init);
    seqan::setDefaultValue(parser, "m",   options.k_max);
    seqan::setDefaultValue(parser, "d",   options.delta_k);

    // Setup hidden options
    setHiddenOptions(parser, true, options);
}


// ==========================================================================
// Function checkInput()
// ==========================================================================

ArgumentParser::ParseResult checkInput(AssemblyOptions & options){

    ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

    if (options.prefix != "." && !exists(options.prefix))
    {
        std::cerr << "ERROR: Path to sample direcotories \'" << options.prefix << "\' does not exist." << std::endl;
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
    std::cerr << "    \033[1massemble\033[0m        Filter, clip and assemble unmapped reads from a sample." << std::endl;
    std::cerr << "    \033[1mmerge\033[0m           Generate supercontigs from a colored compacted de Bruijn Graph." << std::endl;
    std::cerr << "    \033[1mmultik\033[0m          Multi-k framework for a colored compacted de Bruijn Graph." << std::endl;
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
    cout << "temp-path          : " << options.tempPath                 << endl;
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
