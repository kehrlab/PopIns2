/*!
* \file     src/popins2_merge.h
* \brief    Main file for creating a colored compacted de Bruijn Graph out of multiple unitig FASTA files from the
*           popins_single module.
*/
#ifndef POPINS2_MERGE_H_
#define POPINS2_MERGE_H_


#include "argument_parsing.h"           /* seqAn argument parser */
#include "ColoredDeBruijnGraph.h"



/*!
* \fn       int popins2_merge(int argc, char const ** argv)
* \brief    Main function for creating a (merged) colored compacted de Bruijn Graph.
* \return   0 for success, 1 for error
*/
// =========================
// Main
// =========================
int popins2_merge(int argc, char const *argv[]){

    // ==============================
    // Argument Parser
    // ==============================
    CCDBG_Build_opt ccdbg_build_opt;
    CCDBG_Build_opt* p_ccdbg_build_opt = &ccdbg_build_opt;
    MergeOptions mo(p_ccdbg_build_opt);

    seqan::ArgumentParser::ParseResult res = parseCommandLine(mo, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "seqan::ArgumentParser::PARSE_ERROR in module popins2 merge" << endl;
        return 1;
    }


    // ==============================
    // Run graph functions
    // ==============================
    ExtendedCCDBG exg(ccdbg_build_opt.k, ccdbg_build_opt.g);
    cout << "[PROGRESS] Building CCDBG..." << endl;
    exg.buildGraph(ccdbg_build_opt);

    cout << "[PROGRESS] Simplifying CCDBG..." << endl;
    exg.simplify(ccdbg_build_opt.deleteIsolated, ccdbg_build_opt.clipTips, ccdbg_build_opt.verbose);

    cout << "[PROGRESS] ColorMapping CCDBG..." << endl;
    exg.buildColors(ccdbg_build_opt);

    cout << "[PROGRESS] Writing CCDBG..." << endl;
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);


    // TEST START
    /*
    exg.init_ids();
    //exg.connected_components(ccdbg_build_opt);

    for (auto &unitig : exg){
        exg.DFS_Init(unitig, ccdbg_build_opt.verbose);
        exg.DFS_cleaner_seen_only();
    }
    */
    // TEST END




    return 0;
}






#endif /*POPINS2_MERGE_H_*/
