/*!
* \file    src/popins2_single.h
* \brief   Main file for creating a single sample compacted de Bruijn Graph.
*
*/
#ifndef POPINS2_SINGLE_H_
#define POPINS2_SINGLE_H_


#include "CompactedDBG_Graph_extension.h"
#include "argument_parsing.h"           /* seqAn argument parser */



/*!
* \fn       int popins2_single(int argc, char const ** argv)
* \brief    Main function for creating a single sample compacted de Bruijn Graph.
* \return   0 for success, 1 for error
*/
// =========================
// Main
// =========================
int popins2_single(int argc, char const *argv[]){

    // ==============================
    // Argument Parser
    // ==============================
    CDBG_Build_opt cdbg_build_opt;
    CDBG_Build_opt* p_cdbg_build_opt = &cdbg_build_opt;
    SingleOptions so(p_cdbg_build_opt);

    seqan::ArgumentParser::ParseResult res = parseCommandLine(so, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "seqan::ArgumentParser::PARSE_ERROR in module popins2 single" << endl;
        return 1;
    }

    // ==============================
    // Run graph functions
    // ==============================

    ExtendedCDBG xg(cdbg_build_opt.k, cdbg_build_opt.g);
    cout << "[PROGRESS] Building CDBG..." << endl;
    xg.build(cdbg_build_opt);
    cout << "[PROGRESS] Simplifying CDBG..." << endl;
    //cout << xg.size() << endl;
    xg.simplify(cdbg_build_opt.deleteIsolated, cdbg_build_opt.clipTips, cdbg_build_opt.verbose);
    //cout << xg.size() << endl;

    //xg.write("pre_"+cdbg_build_opt.prefixFilenameOut, cdbg_build_opt.nb_threads, cdbg_build_opt.outputGFA, cdbg_build_opt.verbose);

    // TEST start
    xg.init_ids();

    xg.annotate_kmer_coverage(cdbg_build_opt.filename_in);

    /*
    // PRINT kmer cov
    for (auto um : xg){
        DataExtension* de = um.getData();
        cout << de->getID() << ": "; prettyprint::print(de->kmer_cov); cout << endl;
    }
    */

    // xg.write("pre_"+cdbg_build_opt.prefixFilenameOut, cdbg_build_opt.nb_threads, cdbg_build_opt.outputGFA, cdbg_build_opt.verbose);


    xg.small_bubble_removal(cdbg_build_opt.verbose);
    //cout << xg.size() << endl;

    // TEST end

    cout << "[PROGRESS] Writing CDBG..." << endl;
    xg.write(cdbg_build_opt.prefixFilenameOut, cdbg_build_opt.nb_threads, cdbg_build_opt.outputGFA, cdbg_build_opt.verbose);

    //xg.id2csv("post_"+cdbg_build_opt.prefixFilenameOut);    // TODO: make this an cmd option


    return 0;
}







#endif /*POPINS2_SINGLE_H_*/
