/*!
* \file    src/ColoredCDBG_Graph_extension.h
* \brief   Library for a colored compacted de Bruijn Graph using unitig extension.
*
*/
#ifndef COLORED_DE_BRUIJN_GRAPH_
#define COLORED_DE_BRUIJN_GRAPH_


#include <seqan/seq_io.h>
#include <seqan/misc/union_find.h>


#include "UnitigExtension.h"
#include "Traceback.h"
#include "Setcover.h"
#include "PathTracker.h"


// TODO exclude prettyprint for release
#include "../../prettyprint/prettyprint.h"


/*!
* \class        ExtendedCCDBG
* \headerfile   src/ColoredDeBruijnGraph.h
* \brief        Struct to store a colored compacted DBG plus unitig extensions.
*/
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    public:
        ExtendedCCDBG(int kmer_length = 31, int minimizer_length = 23);

        // -------------
        // | Functions |
        // -------------

        void init_ids();
        void print_ids();
        bool is_id_init() const {return id_init_status;}

        void init_entropy();
        bool is_entropy_init() const {return this->entropy_init_status;}


        /**
         *          This function removes every unitig with low entropy from the graph.
         * @param   threshold is the minimum entropy a unitig must have to remain in the graph
         * @return  true if successful
         */
        bool remove_low_entropy(const float threshold);

        /**
         *          Compute the connected component for every node.
         * @ref     seqan/include/seqan/misc/union_find.h
         * @return  true if successful
         */
        bool connected_components(const CCDBG_Build_opt &graph_options);

        /**         Compute the amount of distict connected components in the dBG.
         *          Prior to this function, the connected_components() function needs to be ran successfully.
         * @return  number of connected components
         */
        size_t count_connected_components();

        /**
         *          Main function for the DFS procedure.
         * @param   opt is a wrapper object for the program's input parameter.
         * @return  true if successful
         */
        bool merge(const CCDBG_Build_opt &opt, const int min_kmers, const std::string &outdir);

    private:
        // ----------
        // | Member |
        // ----------

        struct GreaterThan {
            bool operator() (const unsigned &lhs, const unsigned &rhs) const {return lhs>rhs;}
        };

        typedef std::multimap<unsigned, unsigned, GreaterThan> neighborsContainer;

        bool id_init_status;
        bool entropy_init_status;

        seqan::UnionFind<unsigned> UF;

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        // -------------
        // | Functions |
        // -------------

        void DFS_cleaner();
        void DFS_cleaner_seen_only();

        /**
        *           Function to determine the entropy of a UnitigColorMap based on the dinucleotide distribution.
        */
        float entropy(const std::string &sequence);

        /**
         *          Function to determine the traversal direction to go, given a unitig u and the previously visited unitig v.
         *          This function tests the predecessors P of u (due to the orientation of u, v could be predecessor or successor)
         *          and returns GO_FORWARD if v is in P or GO_BACKWARD otherwise.
         * @return  returns a direction of {GO_BACKWARD, GO_FORWARD}
         */
        uint8_t whereToGo(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;

        /**
         *          Function to determine how to reach a previously visited node v with respect to the current unitig u.
         *          Reverses the direction of whereToGo().
         * @details This function can also be used to detemine a neighbor's (NBR) orientation with respect to u if um=NBR and src=u.
         *          NBR  --------->           or    NBR         --------->
         *          u          ----------           u     ----------
         * @return  returns a direction of {GO_BACKWARD, GO_FORWARD}
         */
        uint8_t whereFrom(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;

        /**
         *          Initial stage of the DFS algorithm.
         * @return  Traceback object
         */
        Traceback DFS_Init_bidirectional(const UnitigColorMap<UnitigExtension> &ucm,
                                         const bool verbose,
                                         Setcover<> &sc);

        /**
         *          Recursive stage of the DFS algorithm.
         * @param   src_direction indicates the direction to go (from ucm) to reach the previously visited unitig.
         * @return  Traceback object
         */
        Traceback DFS_Visit_NEW(const UnitigColorMap<UnitigExtension> &ucm,
                                const uint8_t src_direction,
                                Setcover<> &sc,
                                const bool verbose);

        /**
         *          Function to outsource a case analysis for different orientations of the current unitig and its
         *          previously visited unitig.
         * @details See function definition for sketch.
         * @param   neighbor a successor or predecessor of ucm, depending on the traversal direction
         */
        void DFS_case_NEW(const UnitigColorMap<UnitigExtension> &ucm,
                          const UnitigColorMap<UnitigExtension> &neighbor,
                          Traceback &tb,
                          Setcover<> &sc,
                          const bool verbose);

        /**
         *          Funtion to sort neighbors v of a unitig u according to a given criterion f.
         * @param   neighbors is a neighbors object of {BackwardCDBG, ForwardCDBG}
         * @param   container is a data structure D to store and sort the pair (f, v.id) for every v.
         *          D has to be initialized with an appropiate Functor for a descending sort of the pairs.
         */
        template <class TNeighborCDBG>
        void sortNeighbors(const UnitigColorMap<UnitigExtension> &ucm,
                           const TNeighborCDBG &neighbors,
                           neighborsContainer &container) const;

        /**
         *          Funtion to sort (start)nodes v of the dBG according to their kmer length.
         * @param   container is a data structure D to store and sort the pair (len(v), v.id) for every v in the dBG.
         *          D has to be initialized with an appropiate Functor for a descending sort of the pairs.
         */
        void sortStartnodes(neighborsContainer &container) const;

        /**
         *          Function to compute the number of common colors between two unitigs.
         * @details Every color is checked for at least one occurence (in any k-mer) in both unitigs.
         *          We count a color c as common if both unitigs have at least one k-mer with labelled with c.
         * @return  number of common colors
         */
        unsigned check_common_colors(const UnitigColorMap<UnitigExtension> &ucm,
                                     const UnitigColorMap<UnitigExtension> &neighbor) const;

};












#endif /*COLORED_DE_BRUIJN_GRAPH_*/
