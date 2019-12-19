/*!
* @file     src/LECC_Finder.h
* @brief    This file contains the class LECC_Finder.
* @details  The LECC_Finder class is responsibly to identify low entropy connected components (LECC)
*           within a Colored Compacted de Bruijn Graph (CCDBG).
*/
#ifndef LECCFINDER_
#define LECCFINDER_

#include <stdlib.h>
#include <time.h>
#include <unordered_map>
#include "ColoredDeBruijnGraph.h"



class LECC_Finder{

    typedef std::pair<unsigned,unsigned> map_element;
    typedef std::unordered_map<unsigned,unsigned> color_map;

public:
    LECC_Finder(ExtendedCCDBG* exg, const float threshold) :    g_(exg),
                                                                threshold_(threshold) {}

    /**
    *               annotate()
    *   @brief      assigns a LECC identifier to every UnitigColorMap
    *   @details    This function scans every unitig u for it's entropy. If the entropy is
    *               lower than a threshold, a LECC identifier is assigned to u. Then a DFS
    *               is started at u and traverses it's neighborhood V as long as the
    *               entropy remains under the threshold. Every discovered unitig in V gets
    *               assigned the same LECC identifier as u.
    *   @return     unsigned NB_LECCs
    */
    unsigned annotate();

    /**
    *               write()
    *   @brief      This function writes a CSV file with unitigs of LECCs.
    *               Every LECC has most likely a unique color.
    *   @return     true if successful
    */
    bool write(const std::string ofname = "ccdbg.lecc.csv");

private:
    ExtendedCCDBG *g_;

    const float threshold_;

    static char const * const hex_characters;

    bool isGraphInit() const {return g_->is_entropy_init() and g_->is_id_init();}

    /**
    *       annotate_recursion(...)
    *       This function is the recursive part of the DFS in annotate().
    */
    void annotate_recursion(UnitigColorMap<UnitigExtension> &ucm, const unsigned LECC__);

    unsigned create_random_color();

    void u2hex(std::string &hex, unsigned dec, const unsigned length=6);

};









#endif /*LECCFINDER_*/
