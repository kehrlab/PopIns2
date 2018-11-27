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
    
    // DEBUG:
    for (auto file : ccdbg_build_opt.filename_seq_in){
        cout << file << endl;
    }
    
    exg.buildGraph(ccdbg_build_opt);

    cout << "[PROGRESS] Simplifying CCDBG..." << endl;
    exg.simplify(ccdbg_build_opt.deleteIsolated, ccdbg_build_opt.clipTips, ccdbg_build_opt.verbose);

    cout << "[PROGRESS] ColorMapping CCDBG..." << endl;
    exg.buildColors(ccdbg_build_opt);

    cout << "[PROGRESS] Writing CCDBG..." << endl;
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);


    // TEST START
    exg.init_ids();
    if (!exg.is_id_init())
        return -1;

    Traceback tb;
    for (auto &unitig : exg){
        tb = exg.DFS_Init(unitig, ccdbg_build_opt.verbose);
        if (tb.recursive_return_status){
            tb.printIds();
            tb.printSeqs();
            // TODO: concat sequences in Traceback instance
        }
        exg.DFS_cleaner_seen_only();
    }
    // TEST END




    return 0;
}






#endif /*POPINS2_MERGE_H_*/
