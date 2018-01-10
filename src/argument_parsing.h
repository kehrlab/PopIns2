/*!
* \file    src/argument_parsing.h
* \brief   Library for the command line parsing, user interaction and directory navigation.
*
*/
#ifndef ARGUMENT_PARSING_H_
#define ARGUMENT_PARSING_H_


// =========================
// Includes
// =========================
#include <seqan/arg_parse.h>
#include <bifrost/CompactedDBG.hpp>     /* compacted dBG data structure from Bifrost API */


// =========================
// Structs
// =========================
/*!
* \class        OptionsWrapper
* \headerfile   src/argument_parsing.h
* \brief        Struct containing the command line arguments.
* \remark       The parameter are (mostly) a copy of Bifrost's parameter input.
* \ref          https://github.com/pmelsted/bfgraph/blob/master/src/Bifrost.cpp
*/
struct OptionsWrapper{

    // MANDATORY
    std::string indir;
    std::string prefix;
    unsigned unique_kmers;
    unsigned non_unique_kmers;

    // OPTIONAL
    unsigned kmer_length;
    unsigned minimizer_length;
    unsigned threads;

    bool clip_tips;
    bool rm_iso;
    bool verbose;

    std::string bloom_filter_file;

    // default constructor
    OptionsWrapper() :
    verbose(false), clip_tips(false), rm_iso(false)
    {}
};


// =========================
// Functions
// =========================
seqan::ArgumentParser::ParseResult parseCommandLine(OptionsWrapper &options, int argc, char const ** argv);
bool detect_indir_files(OptionsWrapper &options, std::vector<std::string> &sample_fastx_names);
void init_graph_options(OptionsWrapper &options, std::vector<std::string> &sample_fastx_names, CDBG_Build_opt &graph_options);




#endif /*ARGUMENT_PARSING_H_*/