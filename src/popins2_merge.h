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
    // Bifrost graph functions
    // ==============================
    std::ostringstream msg;
    ExtendedCCDBG exg(ccdbg_build_opt.k, ccdbg_build_opt.g);

#ifdef DEBUG
    if (ccdbg_build_opt.verbose) {
        std::cout << "Sequence files:" << std::endl;
        for (auto file : ccdbg_build_opt.filename_seq_in){
            cout << file << endl;
        }
    }
    if (ccdbg_build_opt.verbose) {
        std::cout << "Reference files:" << std::endl;
        for (auto file : ccdbg_build_opt.filename_ref_in){
            cout << file << endl;
        }
    }
#endif // DEBUG

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
    exg.merge(ccdbg_build_opt);


    // ==============================
    // ~~DEBUG~~
    //
    // - prints every unitig as a matrix of color vectors
    // - needs: prettyprint
    // ==============================
    /*
    #include "prettyprint.h"
    size_t nb_colors = exg.getNbColors();

    std::vector<std::vector<bool> > vv;
    for (auto &unitig : exg){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        cout << "Unitig: " << ue->getID() << ", Len: " << unitig.len << endl;

        std::vector<bool> v;
        
        for (size_t kmerIter = 0; kmerIter <= unitig.size - exg.getK(); ++kmerIter){
            const const_UnitigColorMap<UnitigExtension> kmer = unitig.getKmerMapping(kmerIter);
            const UnitigColors* colors = kmer.getData()->getUnitigColors(kmer);

            for (size_t colorIter = 0; colorIter != nb_colors; ++colorIter){
                bool hasColor = colors->contains(kmer, colorIter);
                v.push_back(hasColor);
            }
            vv.push_back(v);
            v.clear();
        }

        prettyprint::print(vv, true);
        vv.clear();
    }
    */





    return 0;
}






#endif /*POPINS2_MERGE_H_*/
