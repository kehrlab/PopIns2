/*!
* \file     src/popins2_merge.h
* \brief    Main file for creating a colored compacted de Bruijn Graph out of multiple unitig FASTA files from the
*           popins_single module.
*/
#ifndef POPINS2_MERGE_H_
#define POPINS2_MERGE_H_


#include "argument_parsing.h"           /* seqAn argument parser */
#include "ColoredDeBruijnGraph.h"
#include "LECC_Finder.h"


typedef std::unordered_map<uint64_t, Kmer> jump_map_t;


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
    MergeOptions mo;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(mo, argc, argv);

    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "[popins2 merge] seqan::ArgumentParser::PARSE_ERROR" << endl;
        return 1;
    }

    CCDBG_Build_opt ccdbg_build_opt;
    setupBifrostOptions(mo, ccdbg_build_opt);

    printMergeOptions(mo);

    std::ostringstream msg;

    // ==============================
    // Bifrost
    // ==============================
    ExtendedCCDBG exg(ccdbg_build_opt.k, ccdbg_build_opt.g);

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
    // PopIns2 merge
    // ==============================
    msg.str("");
    msg << "Assigning ID to every unitig";
    printTimeStatus(msg);
    exg.init_ids();

    msg.str("");
    msg << "Assigning entropy to every unitig";
    printTimeStatus(msg);
    exg.init_entropy();

    ExtendedCCDBG* exg_p = &exg;
    LECC_Finder F(exg_p, 0.7f);

    msg.str("");
    msg << "Computing LECCs";
    printTimeStatus(msg);
    unsigned nb_lecc = F.annotate();

    jump_map_t jump_map;
    jump_map_t* jump_map_ptr;

    msg.str("");
    msg << "Computing jump pairs though LECCs";
    printTimeStatus(msg);
    bool find_jumps_successful = F.find_jumps(jump_map, nb_lecc);

    jump_map_ptr = (find_jumps_successful) ? &jump_map : NULL;

    msg.str("");
    msg << "Connecting jump map with CCDBG";
    printTimeStatus(msg);
    exg.set_jump_map(jump_map_ptr);

    if(mo.write_lecc){
        msg.str("");
        msg << "Writing LECCs";
        printTimeStatus(msg);
        F.write();    
    }

    std::ofstream ofs;
    ofs.open(mo.prefixFilenameOut+".contigs.fasta");

    if(ofs.is_open()){
        msg.str("");
        msg << "Traversing paths in CCDBG";
        printTimeStatus(msg);
        exg.traverse(mo.setcover_min_kmers, ofs, mo.write_setcover, mo.prefixFilenameOut);
    }
    else{
        msg.str("");
        msg << "[popins2_merge()] ERROR opening fasta file";
        return 1;
    }
    ofs.close();

    // ==============================
    // Bifrost
    // ==============================
    msg.str("");
    msg << "Writing CCDBG";
    printTimeStatus(msg);
    exg.write(ccdbg_build_opt.prefixFilenameOut, ccdbg_build_opt.nb_threads, ccdbg_build_opt.verbose);


    return 0;
}






#endif /*POPINS2_MERGE_H_*/
