#include <seqan/arg_parse.h>    /* seqan's lib for argparse */
#include <dirent.h>             /* read folder */
#include "argument_parsing.h"   /* self header */


/*!
* \fn          seqan::ArgumentParser::ParseResult parseCommandLine(OptionsWrapper &options, int argc, char const ** argv)
* \brief       Function handles the input parsing.
*/
seqan::ArgumentParser::ParseResult parseCommandLine(OptionsWrapper &options, int argc, char const ** argv){

    // setup arg parser
    seqan::ArgumentParser parser("popins2");
    addDescription(parser, "Population-scale detection of novel-sequence insertions using de Bruijn Graphs");
#if defined VERSION
    setVersion(parser, VERSION);
#endif
    addUsageLine(parser, "\\--indir DIR \\--prefix STRING \\--unique-kmers INT \\--non-unique-kmers INT [OPTIONS] \\fP ");

    // options
    addOption(parser, seqan::ArgParseOption("i", "indir", "Source directory with FASTQ file(s).", seqan::ArgParseArgument::STRING, "DIRECTORY"));
    setRequired(parser, "indir", true);

    addOption(parser, seqan::ArgParseOption("o", "prefix", "Prefix for the GFA file", seqan::ArgParseArgument::STRING, "TEXT"));
    setRequired(parser, "prefix", true);

    addOption(parser, seqan::ArgParseOption("n", "unique-kmers", "Amount of unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "unique-kmers", true);
    setMinValue(parser, "unique-kmers", "1");

    addOption(parser, seqan::ArgParseOption("N", "non-unique-kmers", "Amount of non-unique kmers.", seqan::ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "non-unique-kmers", true);
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

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Print more output"));

    addOption(parser, seqan::ArgParseOption("c", "clip-tips", "Remove ends of the dBG with length <k"));

    addOption(parser, seqan::ArgParseOption("r", "remove-isolated", "Remove single contigs with length <k"));

    addOption(parser, seqan::ArgParseOption("f", "bloom-filter-file", "Filename for a binary file storing the bloom filter data structure.", seqan::ArgParseArgument::STRING, "TEXT"));

   // parse cmd line
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values to options struct
    getOptionValue(options.indir, parser, "indir");
    getOptionValue(options.prefix, parser, "prefix");
    getOptionValue(options.unique_kmers, parser, "unique-kmers");
    getOptionValue(options.non_unique_kmers, parser, "non-unique-kmers");

    getOptionValue(options.kmer_length, parser, "kmer-length");
    getOptionValue(options.minimizer_length, parser, "minimizer-length");
    getOptionValue(options.threads, parser, "threads");

    options.verbose = isSet(parser, "verbose");
    options.clip_tips = isSet(parser, "clip-tips");
    options.rm_iso = isSet(parser, "remove-isolated");

    getOptionValue(options.bloom_filter_file, parser, "bloom-filter-file");

    // if parse+extraction was successful
    return seqan::ArgumentParser::PARSE_OK;
}


/*!
* \fn           bool init_graph_options(OptionsWrapper &options, CDBG_Build_opt &graph_options)
* \brief        Function parsing the program options to a wrapper struct to initialize a Bifrost graph.
* \details      The general question might arises, why there are two options structs. The global options struct
*               takes care about forcing mandatory inputs and limiting parameters to a certain range. It also creates
*               a help page for the program.
*               The graph options struct wraps all parameters that are exclusively necessary for the graph
*               initialization.
* \return       EXIT_SUCCESS==0 or EXIT_FAILURE==1
*/
bool detect_indir_files(OptionsWrapper &options, std::vector<std::string> &sample_fastx_names){
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(options.indir))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastq = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastq!="." && current_fastq!=".."){
                if (options.verbose)
                    cout << "FILE DETECTED: " << current_fastq << endl;
                sample_fastx_names.push_back(options.indir+current_fastq);
            }
        }
        closedir(dir);
        return EXIT_SUCCESS;
    }
    else {
        // could not open directory
        perror ("");
        return EXIT_FAILURE;
    }
}


/*!
* \fn           bool init_graph_options(OptionsWrapper &options, CDBG_Build_opt &graph_options)
* \brief        Function parsing the program options to a wrapper struct to initialize a Bifrost graph.
* \details      The general question might arises, why there are two options structs. The global options struct
*               takes care about forcing mandatory inputs and limiting parameters to a certain range. It also creates
*               a help page for the program.
*               The graph options struct wraps all parameters that are exclusively necessary for the graph
*               initialization.
*/
void init_graph_options(OptionsWrapper& options, std::vector<std::string> &sample_fastx_names, CDBG_Build_opt& graph_options){
    graph_options.fastx_filename_in = sample_fastx_names;
    graph_options.k = options.kmer_length;
    graph_options.nb_unique_kmers = options.unique_kmers;
    graph_options.nb_non_unique_kmers = options.non_unique_kmers;
    graph_options.prefixFilenameOut = options.prefix;
    graph_options.outFilenameBBF = options.bloom_filter_file;
    graph_options.nb_threads = options.threads;
    graph_options.clipTips = options.clip_tips;
    graph_options.deleteIsolated = options.rm_iso;
}


/*!
* \fn           bool check_ProgramOptions(CDBG_Build_opt &graph_options)
* \brief        Function to check whether the settings are appropriate to construct a CDBG.
* \remark       This function is copied from Bifrost.
* \ref          https://github.com/pmelsted/bfgraph/blob/master/src/Bifrost.cpp
* \return       true if successful
*/
bool check_ProgramOptions(CDBG_Build_opt &graph_options) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    if (graph_options.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
        ret = false;
    }
    else if (graph_options.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << endl;
        ret = false;
    }

    if (graph_options.read_chunksize <= 0) {

        cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (graph_options.k <= 0){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (graph_options.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (graph_options.g <= 0){

        cerr << "Error: Length g of minimizers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (graph_options.g >= MAX_KMER_SIZE){

        cerr << "Error: Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (graph_options.nb_unique_kmers <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!graph_options.reference_mode && (graph_options.nb_non_unique_kmers <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!graph_options.reference_mode && (graph_options.nb_non_unique_kmers > graph_options.nb_unique_kmers)){

        cerr << "Error: The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        ret = false;
    }

    if (graph_options.nb_bits_unique_kmers_bf <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!graph_options.reference_mode && (graph_options.nb_bits_non_unique_kmers_bf <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (graph_options.outFilenameBBF.length() != 0){

        FILE* fp = fopen(graph_options.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "Error: Could not open file for writing output Blocked Bloom filter: " << graph_options.outFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    if (graph_options.inFilenameBBF.length() != 0){

        FILE* fp = fopen(graph_options.inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "Error: Could not read file input Blocked Bloom filter: " << graph_options.inFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    graph_options.filenameOut = graph_options.prefixFilenameOut + (graph_options.outputGFA ? ".gfa" : ".fasta");

    FILE* fp = fopen(graph_options.filenameOut.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << graph_options.filenameOut << endl;
        ret = false;
    }
    else fclose(fp);

    if (graph_options.fastx_filename_in.size() == 0) {

        cerr << "Error: Missing FASTA/FASTQ input files" << endl;
        ret = false;
    }
    else {

        struct stat stFileInfo;
        vector<string>::const_iterator it;
        int intStat;

        for(it = graph_options.fastx_filename_in.begin(); it != graph_options.fastx_filename_in.end(); ++it) {

            intStat = stat(it->c_str(), &stFileInfo);

            if (intStat != 0) {
                cerr << "Error: File not found, " << *it << endl;
                ret = false;
            }
        }
    }

    return ret;
}

