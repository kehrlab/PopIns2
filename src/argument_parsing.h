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
* \fn       inline vector<string> getFilesFromDir(string &path)
* \brief    Function returnes a vector of all files in a given folder (path)
* \remark   Only works for UNIX so far.
*/
inline vector<string> getFilesFromDir(string &path){
    vector<string> sample_fastx_names;      // all file names in --indir with full path
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(path))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastx!="." && current_fastx!="..")
                sample_fastx_names.push_back(path+current_fastx);
        }
        closedir(dir);
    }
    else {
        // could not open directory
        perror ("");
    }

    return sample_fastx_names;
}


/*!
* \fn       inline seqan::ArgumentParser::ParseResult parseCompactedDBGOptions(CDBG_Build_opt &graph_opt, int argc, char const ** argv)
* \brief    Function handles the input parsing for CompactedDBG.
*/
inline seqan::ArgumentParser::ParseResult parseCompactedDBGOptions(CDBG_Build_opt &graph_opt, int argc, char const ** argv){

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


/*!
* \fn       inline seqan::ArgumentParser::ParseResult parseColoredCDBGOptions(CCDBG_Build_opt &graph_opt, int argc, char const ** argv, string &source_path)
* \brief    Function handles the input parsing for ColoredCDBG.
*/
inline seqan::ArgumentParser::ParseResult parseColoredCDBGOptions(CCDBG_Build_opt &graph_opt, int argc, char const ** argv, string &source_path){

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

    addOption(parser, seqan::ArgParseOption("F", "tmp-file-dir", "Source directory with intermediate FASTA file(s) from the popins2 single module.", seqan::ArgParseArgument::STRING, "DIRECTORY"));
    setRequired(parser, "tmp-file-dir", true);

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
    seqan::getOptionValue(source_path, parser, "file-dir");

    // ================================================================================
    // get all FASTA files in tmp directory, that came from CompactedDBG.write()
    // ================================================================================
    string indir;
    seqan::getOptionValue(indir, parser, "tmp-file-dir");
    graph_opt.filename_seq_in = getFilesFromDir(indir);

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