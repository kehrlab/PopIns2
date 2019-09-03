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
    // Bifrost
    // ==============================
    std::ostringstream msg;
    ExtendedCCDBG exg(ccdbg_build_opt.k, ccdbg_build_opt.g);

    // ==============================
    // Bifrost graph functions
    // ==============================
    msg.str("");
    msg << "Building CCDBG";
    printTimeStatus(msg);
    exg.buildGraph(ccdbg_build_opt);

    msg.str("");
    msg << "Simplifying CCDBG";
    printTimeStatus(msg);
    exg.simplify(ccdbg_build_opt.deleteIsolated, ccdbg_build_opt.clipTips, ccdbg_build_opt.verbose);

    msg.str("");
    msg << "ColorMapping CCDBG";
    printTimeStatus(msg);
    exg.buildColors(ccdbg_build_opt);

    // ==============================
    // Bifrost
    // ==============================
    msg.str("");
    msg << "Writing CCDBG";
    printTimeStatus(msg);
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);

    // ==============================
    // Merge specific functions
    // ==============================
    msg.str("");
    msg << "Assigning unitig IDs";
    printTimeStatus(msg);
    exg.init_ids();

    msg.str("");
    msg << "Traversing paths";
    printTimeStatus(msg);
    exg.merge(ccdbg_build_opt, mo.min_kmers);



    return 0;
}






#endif /*POPINS2_MERGE_H_*/
