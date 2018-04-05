/*!
* \file     src/popins2_merge.h
* \brief    Main file for creating a colored compacted de Bruijn Graph out of multiple unitig FASTA files from the
*           popins_single module.
*/
#ifndef POPINS2_MERGE_H_
#define POPINS2_MERGE_H_


#include "argument_parsing.h"           /* seqAn argument parser */
#include "ColoredCDBG_Graph_extension.h"

#include "../../prettyprint/prettyprint.h"      // my debug headder



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
    string source_path;
    seqan::ArgumentParser::ParseResult res = parseColoredCDBGOptions(ccdbg_build_opt, argc, argv, source_path);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "seqan::ArgumentParser::PARSE_ERROR in popins2_merge" << endl;
        return 1;
    }


    // ==============================
    // Run graph functions
    // ==============================
    ExtendedCCDBG exg(ccdbg_build_opt.k, ccdbg_build_opt.g);
    cout << "[PROGRESS] Building CCDBG..." << endl;
    exg.build(ccdbg_build_opt);
    cout << "[PROGRESS] ColorMapping CCDBG..." << endl;
    exg.mapColors(ccdbg_build_opt);

    cout << "[PROGRESS] Writing CCDBG..." << endl;
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);


    // TEST START

    exg.init_ids();
    exg.connected_components(ccdbg_build_opt);

    // WARNING: add a bool whether dfs/bfs (=xfs) variables are on initial state
    for (auto &ucm : exg){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        if (ue->getID() == 3)    // test case: unitig ID 3 as source
            exg.dfs(ucm);
    }

    exg.init_kmer_cov();
    exg.annotate_kmer_coverage(getFilesFromDir(source_path));

    // print k-mer coverage
    /*
    for (auto &ucm : exg){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        cout << "[[" << ue->getID() << "]]";
        prettyprint::print(ue->kmer_coverage);
        cout << endl;
    }
    */

    exg.small_bubble_removal();

    // TEST END







    return 0;
}







#endif /*POPINS2_MERGE_H_*/
