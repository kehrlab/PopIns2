/*!
* \file    src/argument_parsing.h
* \brief   Library for the command line parsing, user interaction and directory navigation.
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
        matepairFile(""), referenceFile(""), prefix("."), sampleID(""),
        kmerLength(47), humanSeqs(maxValue<int>()), threads(1), memory("768M"), use_velvet(false)
    {}
};


struct MergeOptions {

    CCDBG_Build_opt* ccdbg_build_opt;

    string outdir = "";
    signed min_kmers = -1;

    MergeOptions () :
        ccdbg_build_opt(nullptr)
    {}

    MergeOptions (CCDBG_Build_opt* g) :
        ccdbg_build_opt(g)
    {}
};


// =========================
// Option transfer functions
// =========================

bool getOptionValues(AssemblyOptions & options, ArgumentParser const & parser)
{
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
    // Transfer program arguments from command line object to Bifrost's graph options object.
    /* This design recommendation by seqan is absolutely overengineered here, but it leaves the possibility to extend
       the cmd line options object with variables that do not belong to the Bifrost graph. */

    // ---------- Setup merge parameter ----------
    if (seqan::isSet(parser, "outdir")){
        seqan::getOptionValue(options.outdir, parser, "outdir");
        checkPathSyntax(options.outdir);
    }
    else {options.outdir = "";} // default is set here

    if (seqan::isSet(parser, "min-kmers")) seqan::getOptionValue(options.min_kmers, parser, "min-kmers");
    else {options.min_kmers = 64;} // default is set here

    // ---------- Setup graph parameter ----------
    if (seqan::isSet(parser, "verbose")) seqan::getOptionValue(options.ccdbg_build_opt->verbose, parser, "verbose");
  //if (seqan::isSet(parser, "nb_unique_kmers")) seqan::getOptionValue(options.ccdbg_build_opt->nb_unique_kmers, parser, "nb_unique_kmers");
  //if (seqan::isSet(parser, "nb_non_unique_kmers")) seqan::getOptionValue(options.ccdbg_build_opt->nb_non_unique_kmers, parser, "nb_non_unique_kmers");
    if (seqan::isSet(parser, "kmer-length")) seqan::getOptionValue(options.ccdbg_build_opt->k, parser, "kmer-length");
    if (seqan::isSet(parser, "minimizer-length")) seqan::getOptionValue(options.ccdbg_build_opt->g, parser, "minimizer-length");
    if (seqan::isSet(parser, "threads")) seqan::getOptionValue(options.ccdbg_build_opt->nb_threads, parser, "threads");
    if (seqan::isSet(parser, "input-seq-files")) {
        string indir;
        seqan::getOptionValue(indir, parser, "input-seq-files");
        getFastx(options.ccdbg_build_opt->filename_seq_in, indir, options.ccdbg_build_opt->verbose);
    }
    if (seqan::isSet(parser, "input-ref-files")) {
        string indir;
        seqan::getOptionValue(indir, parser, "input-ref-files");
        getFastx(options.ccdbg_build_opt->filename_ref_in, indir, options.ccdbg_build_opt->verbose);
    }
    if (!seqan::isSet(parser, "input-seq-files") && !seqan::isSet(parser, "input-ref-files")){
        std::cout << "[Argparse Error]: At least one of --input-seq-files or --input-ref-files has to be specified." << std::endl;
        return false;
    }
    if (seqan::isSet(parser, "clip-tips")) seqan::getOptionValue(options.ccdbg_build_opt->clipTips, parser, "clip-tips");
    if (seqan::isSet(parser, "del-isolated")) seqan::getOptionValue(options.ccdbg_build_opt->deleteIsolated, parser, "del-isolated");

    // This name defines the prefix for the GFA and bfg_colors file. It is currently not a parameter of the parser.
    if (options.outdir == ""){
        options.ccdbg_build_opt->prefixFilenameOut = "ccdbg";
    }
    else {
        options.ccdbg_build_opt->prefixFilenameOut = options.outdir+"ccdbg";
    }

    return true;
}


// =========================
// Hide options functions
// =========================

void
setHiddenOptions(ArgumentParser & parser, bool hide, AssemblyOptions &)
{
   hideOption(parser, "matePair", hide);
   hideOption(parser, "kmerLength", hide);
}


void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, MergeOptions &){

    seqan::hideOption(parser, "g", hide);
    seqan::hideOption(parser, "m", hide);
  //seqan::hideOption(parser, "n", hide);
  //seqan::hideOption(parser, "N", hide);
}


// ==========================================================================
// Functions setupParser()
// ==========================================================================

void setupParser(ArgumentParser & parser, AssemblyOptions & options)
{
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


/*!
* \fn       void setupParser(seqan::ArgumentParser &parser, MergeOptions &options)
* \brief    Function handles the input parsing for ColoredCDBG.
*/
void setupParser(seqan::ArgumentParser &parser, MergeOptions &options){

    // Setup meta-information
    seqan::setShortDescription(parser, "Build a colored and compacted de Bruijn Graph (CCDBG)");
    seqan::setVersion(parser, VERSION);
    seqan::setDate(parser, DATE);
    seqan::addUsageLine(parser, "\\--input-{seq|ref}-files DIR [OPTIONS] \\fP ");

    // Setup options
    seqan::addSection(parser, "Input/output options");
    seqan::addOption(parser, seqan::ArgParseOption("s", "input-seq-files", "Source directory with FASTA/Q files", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("r", "input-ref-files", "Source directory with reference FASTA/Q files (will not be filtered)", seqan::ArgParseArgument::STRING, "DIR"));
    seqan::addOption(parser, seqan::ArgParseOption("o", "outdir", "Specify a directory for the output files. Default: current working directory.", seqan::ArgParseArgument::STRING, "DIR"));

    seqan::addSection(parser, "Algorithm options");
  //seqan::addOption(parser, seqan::ArgParseOption("n", "unique-kmers", "Amount of unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
  //seqan::addOption(parser, seqan::ArgParseOption("N", "non-unique-kmers", "Amount of non-unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("k", "kmer-length", "K-mer length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("g", "minimizer-length", "Minimizer-length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("v", "verbose", "Print more output"));
    seqan::addOption(parser, seqan::ArgParseOption("i", "clip-tips", "Clip tips shorter than k k-mers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("d", "del-isolated", "Delete isolated contigs shorter than k k-mers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("m", "min-kmers", "Minimum amount of novel k-mers to include a path into the set cover", seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addSection(parser, "Compute resource options");
    seqan::addOption(parser, seqan::ArgParseOption("t", "threads", "Amount of threads for parallel processing", seqan::ArgParseArgument::INTEGER, "INT"));

    // Setup option constraints
    // input (-s/-r) is constrained at getOptionValues()
  //seqan::setRequired(parser, "input-seq-files", true);
  //seqan::setRequired(parser, "output-file", true);

    seqan::setDefaultValue(parser, "k", "31");
    seqan::setDefaultValue(parser, "g", "23");
    seqan::setDefaultValue(parser, "threads", "1");
    seqan::setDefaultValue(parser, "m", "64");     // only for the help print

  //seqan::setMinValue(parser, "unique-kmers", "1");
  //seqan::setMinValue(parser, "non-unique-kmers", "1");
    seqan::setMinValue(parser, "k", "1");
    seqan::setMaxValue(parser, "k", "63");
    seqan::setMinValue(parser, "g", "1");
    seqan::setMaxValue(parser, "g", "62");
    seqan::setMinValue(parser, "t", "1");
    seqan::setMinValue(parser, "m", "1");

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

ArgumentParser::ParseResult checkInput(AssemblyOptions & options)
{
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


ArgumentParser::ParseResult checkInput(MergeOptions & options)
{
    ArgumentParser::ParseResult res = ArgumentParser::PARSE_OK;

    if (options.min_kmers < 0){
        std::cerr << "ERROR: Minimum amount of kmers \'-m/--min-kmers\' can not be less than zero." << std::endl;
        res = ArgumentParser::PARSE_ERROR;
    }

    return res;
}


// =========================
// Parsing functions
// =========================

void printHelp(char const * name)
{
    std::cerr << "Population-scale detection of non-reference sequence insertions using colored de Bruijn Graphs" << std::endl;
    std::cerr << "================================================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " COMMAND\033[0m [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mCOMMAND\033[0m" << std::endl;
    std::cerr << "    \033[1massemble\033[0m        Crop unmapped reads from a bam file and assemble them." << std::endl;
    std::cerr << "    \033[1mmerge\033[0m           Merge many samples into a colored compacted de Bruijn Graph." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << VERSION << ", Date: " << DATE << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}


// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

/*!
* \fn       template<typename TOptions> seqan::ArgumentParser::ParseResult parseCommandLine(TOptions &options, int argc, char const ** argv)
* \brief    Meta-function to setup module based option forwarding
* \remark   Taken from PopIns.
* \return   returns a seqan::ArgumentParser::ParseResult
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