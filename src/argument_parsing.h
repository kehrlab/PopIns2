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

struct MergeOptions {

    CCDBG_Build_opt* ccdbg_build_opt;

    MergeOptions () :
        ccdbg_build_opt(nullptr)
    {}

    MergeOptions (CCDBG_Build_opt* g) :
        ccdbg_build_opt(g)
    {
        // outputColors was automatically set to true during construction time of g
    }
};


// =========================
// Option transfer functions
// =========================

bool getOptionValues(MergeOptions &options, seqan::ArgumentParser &parser){
    // Setup graph object
    if (seqan::isSet(parser, "verbose")) seqan::getOptionValue(options.ccdbg_build_opt->verbose, parser, "verbose");
    if (seqan::isSet(parser, "output-file")) seqan::getOptionValue(options.ccdbg_build_opt->prefixFilenameOut, parser, "output-file");
  //if (seqan::isSet(parser, "nb_unique_kmers")) seqan::getOptionValue(options.ccdbg_build_opt->nb_unique_kmers, parser, "nb_unique_kmers");
  //if (seqan::isSet(parser, "nb_non_unique_kmers")) seqan::getOptionValue(options.ccdbg_build_opt->nb_non_unique_kmers, parser, "nb_non_unique_kmers");
    if (seqan::isSet(parser, "kmer-length")) seqan::getOptionValue(options.ccdbg_build_opt->k, parser, "kmer-length");
    if (seqan::isSet(parser, "minimizer-length")) seqan::getOptionValue(options.ccdbg_build_opt->g, parser, "minimizer-length");
    if (seqan::isSet(parser, "threads")) seqan::getOptionValue(options.ccdbg_build_opt->nb_threads, parser, "threads");
    if (seqan::isSet(parser, "input-seq-files")) {
        string indir;
        seqan::getOptionValue(indir, parser, "input-seq-files");
        // IMPORTANT for debug build: have verbosity check before filename_seq_in assignment
        getFastx(options.ccdbg_build_opt->filename_seq_in, indir, options.ccdbg_build_opt->verbose);
    }
    if (seqan::isSet(parser, "input-ref-files")) {
        string indir;
        seqan::getOptionValue(indir, parser, "input-ref-files");
        // IMPORTANT for debug build: have verbosity check before filename_ref_in assignment
        getFastx(options.ccdbg_build_opt->filename_ref_in, indir, options.ccdbg_build_opt->verbose);
    }
    if (!seqan::isSet(parser, "input-seq-files") && !seqan::isSet(parser, "input-ref-files")){
        std::cout << "[Argparse Error]: At least one of --input-seq-files or --input-ref-files has to be specified." << std::endl;
        return false;
    }
    if (seqan::isSet(parser, "clip-tips")) seqan::getOptionValue(options.ccdbg_build_opt->clipTips, parser, "clip-tips");
    if (seqan::isSet(parser, "del-isolated")) seqan::getOptionValue(options.ccdbg_build_opt->deleteIsolated, parser, "del-isolated");
    if (seqan::isSet(parser, "fasta"))
        options.ccdbg_build_opt->outputGFA = false;

    return true;
}


// =========================
// Hide options functions
// =========================

void setHiddenOptions(seqan::ArgumentParser &parser, bool hide, MergeOptions &){

    seqan::hideOption(parser, "g", hide);
  //seqan::hideOption(parser, "n", hide);
  //seqan::hideOption(parser, "N", hide);
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
  //std::cerr << "    \033[1msingle\033[0m          Build a single sample compacted de Bruijn Graph." << std::endl;
    std::cerr << "    \033[1mmerge\033[0m           Merge many samples into a colored compacted de Bruijn Graph." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << VERSION << ", Date: " << DATE << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
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
    seqan::addUsageLine(parser, "\\--input-{seq|ref}-files DIR \\--output-file STRING [OPTIONS] \\fP ");

    // Setup options
    seqan::addSection(parser, "Input/output options");
    seqan::addOption(parser, seqan::ArgParseOption("s", "input-seq-files", "Source directory with FASTA/Q files", seqan::ArgParseArgument::STRING, "DIRECTORY"));
    seqan::addOption(parser, seqan::ArgParseOption("r", "input-ref-files", "Source directory with reference FASTA/Q files (will not be filtered)", seqan::ArgParseArgument::STRING, "DIRECTORY"));
    seqan::addOption(parser, seqan::ArgParseOption("o", "output-file", "Prefix for the output file", seqan::ArgParseArgument::STRING, "TEXT"));

    seqan::addSection(parser, "Algorithm options");
  //seqan::addOption(parser, seqan::ArgParseOption("n", "unique-kmers", "Amount of unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
  //seqan::addOption(parser, seqan::ArgParseOption("N", "non-unique-kmers", "Amount of non-unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("k", "kmer-length", "K-mer length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("g", "minimizer-length", "Minimizer-length for the dBG construction", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::addOption(parser, seqan::ArgParseOption("v", "verbose", "Print more output"));
    seqan::addOption(parser, seqan::ArgParseOption("i", "clip-tips", "Clip tips shorter than k k-mers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("d", "del-isolated", "Delete isolated contigs shorter than k k-mers in length"));
    seqan::addOption(parser, seqan::ArgParseOption("a", "fasta", "Output file is in FASTA format instead of GFA"));

    seqan::addSection(parser, "Compute resource options");
    seqan::addOption(parser, seqan::ArgParseOption("t", "threads", "Amount of threads for parallel processing", seqan::ArgParseArgument::INTEGER, "INT"));

    // Setup option constraints
    // At least one argument of input-seq-files or input-ref-files has to be specified!
  //seqan::setRequired(parser, "input-seq-files", true);
    seqan::setRequired(parser, "output-file", true);

    seqan::setDefaultValue(parser, "k", "31");
    seqan::setDefaultValue(parser, "g", "23");
    seqan::setDefaultValue(parser, "threads", "1");

  //seqan::setMinValue(parser, "unique-kmers", "1");
  //seqan::setMinValue(parser, "non-unique-kmers", "1");
    seqan::setMinValue(parser, "k", "1");
    seqan::setMaxValue(parser, "k", "63");
    seqan::setMinValue(parser, "g", "1");
    seqan::setMaxValue(parser, "g", "62");
    seqan::setMinValue(parser, "t", "1");

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}


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

    // Transfer program arguments from command line object to Bifrost's graph options object.
    /* This design recommendation by seqan is absolutely overengineered here, but it leaves the possibility to extend
       the cmd line options object with variables that do not belong to the Bifrost graph. */
    if (getOptionValues(options, parser) == false)
        return seqan::ArgumentParser::PARSE_HELP;

    return res;
}






#endif /*ARGUMENT_PARSING_H_*/