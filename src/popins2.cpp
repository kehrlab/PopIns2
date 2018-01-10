#include <iostream>
#include <bifrost/CompactedDBG.hpp>
#include "argument_parsing.h"
using namespace std;



// ==============================
// Function: main()
// ==============================
int main(int argc, char const *argv[]){

    // ==============================
    // Argument Parser
    // ==============================
    OptionsWrapper options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // ==============================
    // Global program variables
    // ==============================
    // all file names in --indir with full path
    vector<string> sample_fastx_names;

    // ==============================
    // Loop for all FastX files
    // ==============================
    bool isReadSuccessful = detect_indir_files(options, sample_fastx_names);

    // ==============================
    // Initialize graph
    // ==============================
    CDBG_Build_opt graph_options;
    if(isReadSuccessful==EXIT_SUCCESS)
        init_graph_options(options, sample_fastx_names, graph_options);



    // ---------- RUN CDBG ----------
    CompactedDBG<> cdbg(graph_options.k, graph_options.g);
    cout << "[PROGRESS] Building CDBG..." << endl;
    cdbg.build(graph_options);

    // ------- SIMPLIFY CDBG -------
    cout << "[PROGRESS] Simplifying CDBG..." << endl;
    cdbg.simplify(graph_options.deleteIsolated, graph_options.clipTips, graph_options.verbose);

    // --------- WRITE CDBG ---------
    cout << "[PROGRESS] Writing GFA..." << endl;
    cdbg.write(graph_options.prefixFilenameOut, graph_options.nb_threads, true, graph_options.verbose);

    cout << "The DBG has " << cdbg.size() << " unitigs.\n" << endl;
}
