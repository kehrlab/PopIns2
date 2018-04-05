/*!
* \file    src/popins2_single.h
* \brief   Main file for creating a single sample compacted de Bruijn Graph.
*
*/
#ifndef POPINS2_SINGLE_H_
#define POPINS2_SINGLE_H_


#include "argument_parsing.h"           /* seqAn argument parser */
#include <bifrost/CompactedDBG.hpp>



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
    seqan::ArgumentParser::ParseResult res = parseGraphOptions(cdbg_build_opt, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "seqan::ArgumentParser::PARSE_ERROR in popins2_single()" << endl;
        return 1;
    }


    // ==============================
    // Run graph functions
    // ==============================
    CompactedDBG<> cdbg(cdbg_build_opt.k, cdbg_build_opt.g);
    cout << "[PROGRESS] Building CDBG..." << endl;
    cdbg.build(cdbg_build_opt);
    cout << "[PROGRESS] Simplifying CDBG..." << endl;
    cdbg.simplify(cdbg_build_opt.deleteIsolated, cdbg_build_opt.clipTips, cdbg_build_opt.verbose);
    cout << "[PROGRESS] Writing CDBG..." << endl;
    cdbg.write(cdbg_build_opt.prefixFilenameOut, cdbg_build_opt.nb_threads, cdbg_build_opt.verbose);

    if (cdbg_build_opt.verbose)
        cout << "[DEBUG] The CDBG has " << cdbg.size() << " unitigs.\n" << endl;




    return 0;
}







#endif /*POPINS2_SINGLE_H_*/
