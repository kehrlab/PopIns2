/*!
* \file    src/ColoredCDBG_Graph_extension.h
* \brief   Library for a colored compacted de Bruijn Graph using unitig extension.
*
*/
#ifndef COLORED_DE_BRUIJN_GRAPH_
#define COLORED_DE_BRUIJN_GRAPH_

#include <seqan/seq_io.h>                   // needed for seqan::log()
#include "UnitigExtension.h"
#include "Traceback.h"


/**
* @class    ExtendedCCDBG
* @brief    This class implemets a colored compacted dBG.
* @param    UnitigExtension is the metadata stored in each UnitigColorMap of the graph
**/
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    struct Greater {
        template<class T>
        bool operator() (T const &lhs, T const &rhs) const {return lhs>rhs;}
    };

    typedef std::multimap<unsigned, unsigned, Greater> ordered_multimap;

    typedef uint8_t direction_t;

public:
    // --------------------
    // | Member functions |
    // --------------------

    /**         Default constructor.
    * @brief    This constructor has to initiate private member variables.
    **/
    ExtendedCCDBG(int kmer_length = 63, int minimizer_length = 23);

    void init_ids();
    void print_ids();
    bool is_id_init() const {return this->id_init_status;}

    void init_entropy();
    bool is_entropy_init() const {return this->entropy_init_status;}

    /**         This function traverses the graph.
    * @return   1 for successful execution
    **/
    uint8_t traverse();


private:
    // --------------------
    // | Member variables |
    // --------------------

    const static direction_t VISIT_SUCCESSOR   = 0x0;
    const static direction_t VISIT_PREDECESSOR = 0x1;

    bool id_init_status;
    bool entropy_init_status;

    // --------------------
    // | Member functions |
    // --------------------

    /**         Depth first search.
                This recursive function contains the main logic for the traversal of the CCDBG.
    * @return   1 for further traversal, 0 for jump back into parent recursion level
    **/
    uint8_t DFS(const UnitigColorMap<UnitigExtension> &ucm, const direction_t direction, Traceback &tb);

    void reset_dfs_states();

    bool is_startnode(const UnitigColorMap<UnitigExtension> &ucm);

    /**         Get a ranking of the neighbors.
    * @brief    This function iterates over all neighbors with respect to the traversal direction. It then applies
    *           the function get_neighbor_overlap() to every neighbor to determine the color-overlap with the current unitig.
    *           It stores the result in an ordered container.
    * @param    omm is a std::multimap<unsigned, unsigned, Greater> to store neighbors in descending order of their color overlap with ucm
    * @param    ucm is the unitig to compare to
    * @param    neighbors is a ForwardCDBG or BackwardCDBG, its unitigs will be compared to ucm
    * @param    direction is the traversal direction
    **/
    template <typename TNeighbors>
    void rank_neighbors(ordered_multimap &omm, const UnitigColorMap<UnitigExtension> &ucm, const TNeighbors &neighbors, const direction_t direction) const;

    /**         Get the color overlap of two neighbor unitig.
    * @brief    The function isolates the kmers that face each other with respect to the unititgs. Then, it retrieves
    *           the color vectors of both kmers,does an AND operation and counts the intersecton.
    * @param    ucm_to_get_head_from is the unitig to get the head kmer from
    * @param    ucm_to_get_tail_from is the unitig to get the tail kmer from
    **/
    unsigned get_neighbor_overlap(const UnitigColorMap<UnitigExtension> &ucm_to_get_head_from, const UnitigColorMap<UnitigExtension> &ucm_to_get_tail_from) const;

    /**         Computes the entropy for a given string.
     * @brief   If all dimers are equaly distributed, the entropy is high (highly chaotic system),
     *          if certain dimers are prevalent, the entropy is low (highly ordered system).
     * @return  The entropy [0,1] of bi-nucleotides
     */
    float entropy(const std::string &sequence);

    /* DEBUG functions */
    inline void print_unitig_id(const UnitigColorMap<UnitigExtension> &ucm){DataAccessor<UnitigExtension>* da = ucm.getData(); UnitigExtension* data = da->getData(ucm); std::cout << data->getID() << std::endl;}
    inline unsigned get_unitig_id(const UnitigColorMap<UnitigExtension> &ucm){DataAccessor<UnitigExtension>* da = ucm.getData(); UnitigExtension* data = da->getData(ucm); return data->getID();}
    inline void print_unitig_seq(const UnitigColorMap<UnitigExtension> &ucm){std::cout << ucm.mappedSequenceToString() << std::endl;}
    inline std::string get_unitig_seq(const UnitigColorMap<UnitigExtension> &ucm){return ucm.mappedSequenceToString();}
};












#endif /*COLORED_DE_BRUIJN_GRAPH_*/
