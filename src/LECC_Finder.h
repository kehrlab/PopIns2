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

    typedef uint8_t direction_t;

    typedef std::pair<unsigned,unsigned> color_map_element;     // only used for writing a LECC color file
    typedef std::unordered_map<unsigned,unsigned> color_map;    // only used for writing a LECC color file

    typedef std::unordered_map<Kmer, bool, KmerHash> border_map_t;

    typedef std::unordered_map<uint64_t, Kmer> jump_map_t;

public:
    // --------------------
    // | Member functions |
    // --------------------

    /**
    *               Constructor
    *   @param      exg is a pointer to a CCDBG
    *   @param      threshold is the upper bound in [0.25,1] to mask unitigs as low entropy
    */
    LECC_Finder(ExtendedCCDBG* exg, const float threshold) :    g_(exg),
                                                                threshold_(threshold),
                                                                lecc_init_status_(false) {}

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
    *               Every LECC has (most likely) a unique color.
    *   @return     true if successful
    */
    bool write(const std::string ofname = "ccdbg.lecc.csv");


    /**
    *               find_jumps()
    *   @brief      determines the best jumps across a LECC
    *   @detail     For every border kmer this function searches for the best color-matching
    *               partner that is accessible across the LECC.
    *   @return     true if successful
    */
    bool find_jumps(jump_map_t &jump_map, const unsigned nb_leccs);


private:
    // --------------------
    // | Member variables |
    // --------------------

    const static direction_t VISIT_SUCCESSOR   = 0x0;
    const static direction_t VISIT_PREDECESSOR = 0x1;

    ExtendedCCDBG *g_;

    const float threshold_;

    bool lecc_init_status_;

    static char const * const hex_characters;

    // --------------------
    // | Member functions |
    // --------------------

    bool isGraphInit() const {return g_->is_entropy_init() and g_->is_id_init();}


    /**
    *       annotate_recursion()
    *       This function is the recursive part of the bidirectional DFS in annotate().
    *       @param  ucm is the current node
    *       @param  LECC__ is the counter of overall LECCs
    *       @param  d is the direction of the traversal
    */
    void annotate_recursion(UnitigColorMap<UnitigExtension> &ucm, const unsigned LECC__);


    unsigned create_random_color();


    void u2hex(std::string &hex, unsigned dec, const unsigned length=6);


    /**
    *       get_borders()
    *       This function detects the bordering unitigs of a LECC
    *       @param  lecc_id is a LECC identifier to compare to
    *       @param  border_kmers is a container to store the bordering Kmers in (the unitig's kmer facing towards the LECC)
    *       @return true if LECC has borders, false otherwise (happens if LECC is a singleton)
    */
    bool get_borders(border_map_t &border_kmers, const unsigned lecc_id) const;


    /**
    *       color_overlap()
    *       This function calculates the color overlap of two kmers
    *       @param  kmer1 is an instance of the Kmer class containing exactly one kmer
    *       @param  kmer2 is an instance of the Kmer class containing exactly one kmer
    *       @return jaccard index of the color overlap
    */
    float color_overlap(const Kmer &kmer1, const Kmer &kmer2) const;


    void reset_dfs_states() const;


    /**
    *       DFS_init()
    *       Initiates a directed depth first serach trying to find a jump over a LECC to one of the n-best color matches.
    *       @param  ucm is the current state (unitig) of the DFS
    *       @param  border_kmers is a container to store the bordering Kmers in (the unitig's kmer facing towards the LECC).
    */
    void check_accessibility(border_map_t &border_kmers, const Kmer &kmer, const unsigned lecc_id) const;


    /**
    *       DFS()
    *       Directed depth first serach trying to find a jump over a LECC to one of the n-best color matches.
    *       @param  ucm is the current state (unitig) of the DFS
    *       @param  d is the traversal direction
    *       @param  border_kmers is a container to store the bordering Kmers in (the unitig's kmer facing towards the LECC).
    *       @return bool; 0 if error occurred; 1 if everything was alright
    */
    bool DFS(border_map_t &border_kmers, const UnitigColorMap<UnitigExtension> &ucm, const direction_t d) const;

};









#endif /*LECCFINDER_*/
