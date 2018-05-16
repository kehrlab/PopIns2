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
 * \fn      template<typename TGraph> string dye3samples(const UnitigColorMap<UnitigExtension> &ucm, TGraph &g)
 * \brief   This function is used to dye a three samples colored de Bruijn Graph, i.e. write a color encoding
 *          out to the GFA file.
 * \return  string; color value in hex, e.g. "#be0204"
 */
string dye3samples(const UnitigColorMap<UnitigExtension> &ucm, const ExtendedCCDBG &g){

    size_t nb_colors = g.getNbColors();

    if (nb_colors != 3) return string();

    size_t k = g.getK();

    std::bitset<3> flip1stBit(std::string("100"));
    std::bitset<3> flip2ndBit(std::string("010"));
    std::bitset<3> flip3rdBit(std::string("001"));

    std::bitset<3> gelb  (std::string("100"));
    std::bitset<3> rot   (std::string("010"));
    std::bitset<3> blau  (std::string("001"));
    std::bitset<3> orange(std::string("110"));
    std::bitset<3> gruen (std::string("101"));
    std::bitset<3> lila  (std::string("011"));
    std::bitset<3> braun (std::string("111"));

    // For each k-mer in the Unitig
    for (size_t i = 0; i <= ucm.size - k; ++i){

        const const_UnitigColorMap<UnitigExtension> ucm_mapping = ucm.getKmerMapping(i);        // Get the mapping for this k-mer
        const UnitigColors* uc = ucm_mapping.getData()->getUnitigColors(ucm_mapping);           // Get the color set associated with this unitig

        std::bitset<3> colvec(std::string("000"));

        // modify color vector according to all three color sets
        if(uc->contains(ucm_mapping, 0)) colvec |= flip1stBit;
        if(uc->contains(ucm_mapping, 1)) colvec |= flip2ndBit;
        if(uc->contains(ucm_mapping, 2)) colvec |= flip3rdBit;

        // return a color depending on the bitmask of colvec
        if      (colvec == gelb)   return "\tCL:z:#f4ed1f\tC2:z:#f4ed1f";
        else if (colvec == rot )   return "\tCL:z:#f42018\tC2:z:#f42018";
        else if (colvec == blau)   return "\tCL:z:#1401e2\tC2:z:#1401e2";
        else if (colvec == orange) return "\tCL:z:#e69421\tC2:z:#e69421";
        else if (colvec == gruen)  return "\tCL:z:#26c814\tC2:z:#26c814";
        else if (colvec == lila)   return "\tCL:z:#bc2da9\tC2:z:#bc2da9";
        else                       return "\tCL:z:#7c6241\tC2:z:#7c6241";   // braun
    }
}

/*
struct Dye3samples {

    Dye3samples(const ExtendedCCDBG &g) : g_(g) {}

    string operator()()


    private:

        ExtendedCCDBG g_;
};
*/

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
    exg.build(ccdbg_build_opt);
    cout << "[PROGRESS] ColorMapping CCDBG..." << endl;
    exg.mapColors(ccdbg_build_opt);

    cout << "[PROGRESS] Writing CCDBG..." << endl;
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);


    // TEST START

    

    //exg.init_ids();
    //exg.connected_components(ccdbg_build_opt);

    /*
    // WARNING: add a bool whether dfs/bfs (=xfs) variables are on initial state
    for (auto &ucm : exg){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        if (ue->getID() == 3)    // test case: unitig ID 3 as source
            exg.dfs(ucm);
    }
    */

    //exg.init_kmer_cov();
    //exg.annotate_kmer_coverage(getFilesFromDir(mo.source_path)); // NOTE: deleted source_path

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

    //exg.small_bubble_removal();   // NOTE: if needed then needs to be reworked, VERY OLD

    // TEST END







    return 0;
}







#endif /*POPINS2_MERGE_H_*/
