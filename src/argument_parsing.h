/*!
* \file    src/argument_parsing.h
* \brief   Library for the command line parsing, user interaction and directory navigation.
*
*/
#ifndef ARGUMENT_PARSING_H_
#define ARGUMENT_PARSING_H_


#include <bifrost/ColoredCDBG.hpp>

#include <vector>
#include <dirent.h>             /* read folder */
#include <seqan/arg_parse.h>
using namespace std;



// =========================
// Functions
// =========================

/*!
* \fn       template <typename TType> seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier(TType &graph_opt, seqan::ArgumentParser &parser)
* \brief    This function should never be called. Only exists to satisfy the compiler.
*/
template <typename TType>
inline seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier(TType &graph_opt, seqan::ArgumentParser &parser){
    // catch case for corrupted template deduction
    cerr << "ERROR: Undefined template parameter (should be a Bifrost graph type)." << endl;

    return seqan::ArgumentParser::PARSE_ERROR;
}


/*!
* \fn       template <> seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier<CCDBG_Build_opt>(CCDBG_Build_opt &graph_opt, seqan::ArgumentParser &parser)
* \brief    Parse options specific to ColoredCDBG_Build_opt
*/
template <>
inline seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier<CCDBG_Build_opt>(CCDBG_Build_opt &graph_opt, seqan::ArgumentParser &parser){

    // ================================================================================
    // get all FASTA files in input directory, that came from CompactedDBG.write()
    // ================================================================================
    string indir;
    seqan::getOptionValue(indir, parser, "file-dir");
    vector<string> sample_fastx_names;      // all file names in --indir with full path
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(indir))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastx!="." && current_fastx!=".."){
                if (graph_opt.verbose)
                    cout << "FILE DETECTED: " << current_fastx << endl;
                sample_fastx_names.push_back(indir+current_fastx);
            }
        }
        closedir(dir);
    }
    else {
        // could not open directory
        perror ("");
    }
    graph_opt.filename_seq_in = sample_fastx_names;


    return seqan::ArgumentParser::PARSE_OK;
}


/*!
* \fn       template <> seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier<CDBG_Build_opt>(CDBG_Build_opt &graph_opt, seqan::ArgumentParser &parser)
* \brief    Parse options specific to CompactedDBG_Build_opt
*/
template <>
inline seqan::ArgumentParser::ParseResult parseGraphOptionsSpecifier<CDBG_Build_opt>(CDBG_Build_opt &graph_opt, seqan::ArgumentParser &parser){

    // ================================================================================
    // get all FASTA files in input directory, that came from CompactedDBG.write()
    // ================================================================================
    string indir;
    seqan::getOptionValue(indir, parser, "file-dir");
    vector<string> sample_fastx_names;      // all file names in --indir with full path
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(indir))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastx!="." && current_fastx!=".."){
                if (graph_opt.verbose)
                    cout << "FILE DETECTED: " << current_fastx << endl;
                sample_fastx_names.push_back(indir+current_fastx);
            }
        }
        closedir(dir);
    }
    else {
        // could not open directory
        perror ("");
    }
    graph_opt.filename_in = sample_fastx_names;


    return seqan::ArgumentParser::PARSE_OK;
}


/*!
* \fn       template <typename TType> seqan::ArgumentParser::ParseResult parseCommandLine(TType &graph_opt, int argc, char const ** argv)
* \var      T
*           template parameter for the type of graph
* \brief    Function handles the input parsing.
*/
template <typename TType>
inline seqan::ArgumentParser::ParseResult parseGraphOptions(TType &graph_opt, int argc, char const ** argv){

    // ==================
    // setup arg parser
    // ==================

    seqan::ArgumentParser parser("popins2");
    addDescription(parser, "Population-scale detection of novel-sequence insertions using colored de Bruijn Graphs");
#if defined VERSION
    setVersion(parser, VERSION);
#endif
    addUsageLine(parser, "\\--file-dir DIR \\--output-file STRING [OPTIONS] \\fP ");

    // options
    addOption(parser, seqan::ArgParseOption("f", "file-dir", "Source directory with input FASTX file(s).", seqan::ArgParseArgument::STRING, "DIRECTORY"));
    setRequired(parser, "file-dir", true);

    addOption(parser, seqan::ArgParseOption("o", "output-file", "Prefix for the GFA file", seqan::ArgParseArgument::STRING, "TEXT"));
    setRequired(parser, "output-file", true);

    addOption(parser, seqan::ArgParseOption("n", "unique-kmers", "Amount of unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "unique-kmers", "1");

    addOption(parser, seqan::ArgParseOption("N", "non-unique-kmers", "Amount of non-unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "non-unique-kmers", "1");

    addOption(parser, seqan::ArgParseOption("k", "kmer-length", "K-mer length for the dBG construction.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "1");
    setMaxValue(parser, "k", "63");

    addOption(parser, seqan::ArgParseOption("g", "minimizer-length", "Minimizer-length for the dBG construction.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "g", "23");
    setMinValue(parser, "g", "1");
    setMaxValue(parser, "g", "62");

    addOption(parser, seqan::ArgParseOption("t", "threads", "Amount of threads for parallel processing.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", "1");
    setMinValue(parser, "t", "1");

    addOption(parser, seqan::ArgParseOption("i", "clip-tips", "Remove branching ends (tips) of the dBG shorter than k k-mers in length"));

    addOption(parser, seqan::ArgParseOption("d", "del-isolated", "Remove single contigs shorter than k k-mers in length"));

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Print more output"));

   // parse cmd line
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // =====================
    // setup graph options
    // =====================

    // Extract generic option values to graph struct
    graph_opt.verbose = isSet(parser, "verbose");
    graph_opt.clipTips = isSet(parser, "clip-tips");
    graph_opt.deleteIsolated = isSet(parser, "del-isolated");
    seqan::getOptionValue(graph_opt.prefixFilenameOut, parser, "output-file");      // <-- mandatory
    seqan::getOptionValue(graph_opt.nb_unique_kmers, parser, "unique-kmers");
    seqan::getOptionValue(graph_opt.nb_non_unique_kmers, parser, "non-unique-kmers");
    seqan::getOptionValue(graph_opt.k, parser, "kmer-length");
    seqan::getOptionValue(graph_opt.g, parser, "minimizer-length");
    seqan::getOptionValue(graph_opt.nb_threads, parser, "threads");

    // Extract template type specific option values to graph struct
    if (parseGraphOptionsSpecifier(graph_opt, parser) != seqan::ArgumentParser::PARSE_OK)
        return seqan::ArgumentParser::PARSE_ERROR;

    // ==================================
    // catch some common mistakes
    // ==================================

    size_t max_threads = std::thread::hardware_concurrency();
    if (graph_opt.nb_threads > max_threads){
        cerr << "Error: Number of threads has to be smaller than " << max_threads << endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }


    return seqan::ArgumentParser::PARSE_OK;
}




#endif /*ARGUMENT_PARSING_H_*/