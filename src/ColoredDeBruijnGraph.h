/**
 * @file    src/ColoredCDBG_Graph_extension.h
 * @brief   Library for a colored compacted de Bruijn Graph using unitig extension.
 */
#ifndef COLORED_DE_BRUIJN_GRAPH_
#define COLORED_DE_BRUIJN_GRAPH_

#include <seqan/seq_io.h>                   // needed for seqan::log()
#include "UnitigExtension.h"
#include "Traceback.h"
#include "Setcover.h"

#include "debug_macros.h"
#include "prettyprint.h"        // TODO: delete at release


/**
 * @class   ExtendedCCDBG
 * @brief   This class implemets a colored compacted dBG.
 * @param   UnitigExtension is the metadata stored in each UnitigColorMap of the graph
 */
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    struct Greater {
        template<class T>
        bool operator() (T const &lhs, T const &rhs) const {return lhs>rhs;}
    };

    typedef std::multimap<float, unsigned, Greater> ordered_multimap_t;

    typedef std::unordered_map<uint64_t, Kmer> jump_map_t;

    typedef uint8_t direction_t;

public:
    // --------------------
    // | Member functions |
    // --------------------

    /**
     *          Default constructor.
     * @brief   This constructor has to initiate private member variables.
     */
    ExtendedCCDBG(int kmer_length = 63, int minimizer_length = 23) :
        ColoredCDBG<UnitigExtension> (kmer_length, minimizer_length),
        id_init_status(false),
        entropy_init_status(false)
    {}

    void init_ids();
    void print_ids();
    bool is_id_init() const {return this->id_init_status;}

    void init_entropy();
    bool is_entropy_init() const {return this->entropy_init_status;}


    /**
     *          This function traverses the graph.
     * @param   setcover_threshold is the required minimum amount of unseen kmers
     *          to include a path into the final solution.
     * @param   ofs is an outfile stream to write in the contig
     * @param   write_setcover is a boolean indicating whether the IDs of the setcover
     *          should be written to a file.
     * @param   prefixFilenameOut is the prefix for all filenames
     *          (currently only used for Setcover::write)
     * @return  1 successful execution
     *          0 sanity check(s) failed
     */
    uint8_t traverse(const int setcover_threshold, ofstream &ofs, const bool write_setcover, const string prefixFilenameOut);


    void set_jump_map(jump_map_t* m) {_jump_map_ptr = m;}


private:
    // --------------------
    // | Member variables |
    // --------------------

    const static direction_t VISIT_SUCCESSOR   = 0x0;
    const static direction_t VISIT_PREDECESSOR = 0x1;

    bool id_init_status;
    bool entropy_init_status;

    jump_map_t *_jump_map_ptr = NULL;

    // --------------------
    // | Member functions |
    // --------------------

    /**
     *          This Depth First Search function contains the main logic
     *          for the traversal of the CCDBG.
     * @param   ucm is the current unitig of the traversal
     * @param   direction is the traversal direction (VISIT_PREDECESSOR/ VISIT_SUCCESSOR)
     * @param   tb is a Traceback instance that concatenates unititg sequences
     * @param   sc is a Setcover instance that keeps track of the overall written unitigs
     * @param   jumped is an indicator if the latest traversal step (to ucm) was
     *          a jumep over a LECC (default: false)
     * @return  1 for further traversal, 0 for jump back into parent recursion level
     */
    uint8_t DFS(const UnitigColorMap<UnitigExtension> &ucm, const direction_t direction, Traceback &tb, Setcover &sc, const bool jumped = false);


    void reset_dfs_states();


    bool is_startnode(const UnitigColorMap<UnitigExtension> &ucm) const;


    /**         Get a ranking of the neighbors.
     * @brief   This function iterates over all neighbors with respect to the
     *          traversal direction. It then applies the function get_neighbor_overlap()
     *          to every neighbor to determine the color-overlap with the current
     *          unitig. It stores the result in an ordered container.
     * @param   omm is a std::multimap<unsigned, unsigned, Greater> to store neighbors
     *          in descending order of their color overlap with ucm
     * @param   ucm is the unitig to compare to
     * @param   neighbors is a ForwardCDBG or BackwardCDBG, its unitigs will be compared to ucm
     * @param   direction is the traversal direction
     */
    template <typename TNeighbors>
    void rank_neighbors(ordered_multimap_t &omm, const UnitigColorMap<UnitigExtension> &ucm, const TNeighbors &neighbors, const direction_t direction);


    /**         Get the color overlap of two neighbor unitig.
     * @brief   The function isolates the kmers that face each other with respect
     *          to the unititgs. Then, it retrieves the color vectors of both
     *          kmers,does an AND operation and counts the intersecton.
     * @param   extract_head is the unitig to get the head kmer from
     * @param   extract_tail is the unitig to get the tail kmer from
     * @return  jaccard index of the colors of the overlapping kmers
     */
    float get_neighbor_overlap(const UnitigColorMap<UnitigExtension> &extract_head, const UnitigColorMap<UnitigExtension> &extract_tail) ;


    /**         Computes the entropy for a given string.
     * @brief   If all dimers are equaly distributed, the entropy is high (highly chaotic system),
     *          if certain dimers are prevalent, the entropy is low (highly ordered system).
     * @return  The entropy [0,1] of bi-nucleotides
     */
    float entropy(const std::string &sequence);


    /**         Determines the direction to continue after a jump over a LECC
     * @brief   After the traversal jumped over a LECC, this function determines in which
     *          direction (w.r.t. the current unitig) the traversal has to continue.
     *          NOTE: This function *might* be redundant since the traversal could just continue
     *          in the same direction as before the jump. But I am not sure if e.g. a palindromic
     *          self-loop can inverse the traversal direction. This function introduces only a
     *          neglectable computational overhead though.
     * @param   partner is the LECC border Kmer that the traversal jumped to
     * @return  a traversal direction (VISIT_SUCCESSOR or VISIT_PREDECESSOR) or 0x2 (ERROR STATE)
     */
    uint8_t post_jump_continue_direction(const UnitigColorMap<UnitigExtension> &ucm) const;


    // -------------------
    // | Debug functions |
    // -------------------

    inline unsigned get_unitig_id(const UnitigColorMap<UnitigExtension> &ucm) const{
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* data = da->getData(ucm);
        return data->getID();
    }


    inline std::string get_unitig_seq(const UnitigColorMap<UnitigExtension> &ucm) const{
        return ucm.mappedSequenceToString();
    }


    inline unsigned get_unitig_lecc(const UnitigColorMap<UnitigExtension> &ucm) const{
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* data = da->getData(ucm);
        return data->getLECC();
    }
};












#endif /*COLORED_DE_BRUIJN_GRAPH_*/
