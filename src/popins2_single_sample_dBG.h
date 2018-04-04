/*!
* \file    src/popins2_single_sample_dBG.h
* \brief   Main file for creating a single sample compacted de Bruijn Graph.
*
*/
#ifndef POPINS2_SINGLE_SAMPLE_DBG_H_
#define POPINS2_SINGLE_SAMPLE_DBG_H_

#include "argument_parsing.h"           /* seqAn argument parser */



/*!
* \fn       int popins_single_sample_dBG(int argc, char const ** argv)
* \brief    Main function for creating a single sample compacted de Bruijn Graph.
* \return   0 for success, 1 for error
*/
// =========================
// Main
// =========================
int popins_single_sample_dBG(int argc, char const *argv[]){
    // ==============================
    // Argument Parser
    // ==============================
    CDBG_Build_opt cdbg_build_opt;
    seqan::ArgumentParser::ParseResult res = parseGraphOptions(cdbg_build_opt, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK){
        cerr << "seqan::ArgumentParser::PARSE_ERROR in popins_single_sample_dBG()" << endl;
        return 1;
    }




    return 0;
}







#endif /*POPINS2_SINGLE_SAMPLE_DBG_H_*/