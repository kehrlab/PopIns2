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
#include "setcover.h"


// TODO exclude prettyprint for release
#include "../../prettyprint/prettyprint.h"


typedef std::vector<std::vector<std::string> >VVSequences;
typedef std::vector<std::string> VSequences;


// =========================
// Structs
// =========================

/*!
* \class        Traceback
* \headerfile   src/ColoredDeBruijnGraph.h
* \brief        Struct to manage the metadata for the DFS traceback.
*/
class Traceback{

private:
    void cutconcat(string &s, const VSequences &path, const size_t k) const;

    VVSequences pathseqs;

public:
    using const_iterator = VVSequences::const_iterator;
    using iterator = VVSequences::iterator;
    using const_reverse_iterator = VVSequences::const_reverse_iterator;

    bool recursive_return_status = false;

    static uint8_t recursion_priority_counter;

    bool write(ofstream &ofs, const size_t k, size_t &counter) const;

    void join(const Traceback &t);

    void rearrange(const Traceback &bw, const Traceback &fw);

    bool empty() const {return pathseqs.empty();}
    void push_back(const VSequences &ps) {pathseqs.push_back(ps);}
    const_iterator cbegin() const { return pathseqs.cbegin(); }
    const_iterator cend() const { return pathseqs.cend(); }
    iterator begin() { return pathseqs.begin(); }
    iterator end() { return pathseqs.end(); }
    const_reverse_iterator crbegin() const {return pathseqs.crbegin();}
    const_reverse_iterator crend() const {return pathseqs.crend();}
};


/*!
* \class        ExtendedCCDBG
* \headerfile   src/ColoredDeBruijnGraph.h
* \brief        Struct to store a colored compacted DBG plus unitig extensions.
*/
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    public:
        ExtendedCCDBG(int kmer_length = 31, int minimizer_length = 23);      // hidden inits! (see definition)

        void init_ids();
        void print_ids();
        bool is_id_init() const {return id_init_status;}

        bool connected_components(const CCDBG_Build_opt &graph_options);
        size_t count_connected_components();
        seqan::UnionFind<unsigned> getUF() const {return UF;}

        bool merge(const CCDBG_Build_opt &opt, const unsigned max_paths);

    private:
        struct GreaterThan {
            bool operator() (const char &lhs, const char &rhs) const {return lhs>rhs;}
        };

        typedef std::multimap<float, unsigned, GreaterThan> neighborsContainer;

        typedef std::pair <float, size_t> keyPair;

        struct keyPair_GreaterThan {
            bool operator() (const keyPair &lhs, const keyPair &rhs) const {
                if (lhs.first != rhs.first) {return lhs.first > rhs.first;}
                else {return lhs.second > rhs.second;}
            }
        };

        typedef std::multimap<keyPair, unsigned, keyPair_GreaterThan> neighborsContainerWithStoredLength;

        bool id_init_status;

        seqan::UnionFind<unsigned> UF;

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        float entropy(const std::string &sequence);

        uint8_t whereToGo(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;
        uint8_t whereFrom(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;

        Traceback DFS_Init(const UnitigColorMap<UnitigExtension> &ucm,
                           Setcover<> &sc,
                           const bool verbose,
                           const unsigned max_paths);

        Traceback DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
                            std::vector<bool> &start_vec,
                            const uint8_t src_direction,
                            Setcover<> &sc,
                            const bool verbose,
                            const unsigned max_paths,
                            const bool started_at_internal);

        void DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
                      const UnitigColorMap<UnitigExtension> &neighbor,
                      std::vector<bool> &start_vec,
                      Traceback &tb,
                      Setcover<> &sc,
                      const bool verbose,
                      const unsigned max_paths,
                      const bool started_at_internal);

        void DFS_cleaner();
        void DFS_cleaner_seen_only();

        bool yStemCheck(const UnitigColorMap<UnitigExtension> &startnode,
                             const bool verbose);

        bool haveCommonColor(const UnitigColorMap<UnitigExtension> &start_ucm,
                                 const UnitigColorMap<UnitigExtension> &ucm,
                                 const uint8_t start_direction,
                                 const bool verbose) const;

        float equalColorbitsRate(const std::vector<bool> &v,
                                 const UnitigColorMap<UnitigExtension> &neighbor) const;

        template <class TContainer> void getSourceNodes(TContainer &m) const;

        void update_start_vec(std::vector<bool> &start_vec,
                              const UnitigColorMap<UnitigExtension> &ucm) const;
        bool is_empty_start_vec(const std::vector<bool> &start_vec) const;


        template <class TNeighborCDBG>
        void sortNeighbors(const TNeighborCDBG &neighbors,
                           const std::vector<bool> &start_vec,
                           neighborsContainer &container) const;

        template <class TNeighborCDBG>
        void sortNeighbors(const TNeighborCDBG &neighbors,
                           const std::vector<bool> &start_vec,
                           neighborsContainerWithStoredLength &container) const;


};












#endif /*COLORED_DE_BRUIJN_GRAPH_*/