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

// TODO exclude prettyprint for release
#include "../../prettyprint/prettyprint.h"

#ifdef DEBUG
typedef std::vector<std::vector<unsigned> > PathSet;
typedef std::vector<unsigned> Path;
#endif // DEBUG

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

    bool recursive_return_status = false;

#ifdef DEBUG
    PathSet ids;
    std::vector<std::vector<bool> > oris;
    std::vector<std::vector<std::string> > seqs;

    void printIds() const;
    void printOris() const;
    void printSeqs() const;
    void printPathSeqs() const;
#endif // DEBUG

    bool write(ofstream &ofs, const size_t k, size_t &counter) const;

    void join(const Traceback &t);

    bool empty() const {return pathseqs.empty();}
    void push_back(const VSequences &ps) {pathseqs.push_back(ps);}
    const_iterator cbegin() const { return pathseqs.cbegin(); }
    const_iterator cend() const { return pathseqs.cend(); }
    iterator begin() { return pathseqs.begin(); }
    iterator end() { return pathseqs.end(); }
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

        bool merge(const CCDBG_Build_opt &opt);

    private:

        bool id_init_status;

        seqan::UnionFind<unsigned> UF;

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        float entropy(const std::string &sequence);

        uint8_t whereToGo(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;
        uint8_t whereFrom(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;

        Traceback DFS_Init(const UnitigColorMap<UnitigExtension> &ucm, const bool verbose);

        Traceback DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
                            const UnitigColorMap<UnitigExtension> &start_ucm,
                            const uint8_t src_direction,
                            const uint8_t start_direction,
                            const bool verbose);

        void DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
                      const UnitigColorMap<UnitigExtension> &neighbor,
                      const UnitigColorMap<UnitigExtension> &start_ucm,
                      Traceback &tb,
                      const uint8_t start_direction,
                      const bool verbose);
        
        void DFS_cleaner();
        void DFS_cleaner_seen_only();

        bool haveCommonColor(const UnitigColorMap<UnitigExtension> &start_ucm,
                                 const UnitigColorMap<UnitigExtension> &ucm,
                                 const uint8_t start_direction, 
                                 const bool verbose) const;

        template <class TContainer> void getSourceNodes(TContainer &m) const;

};












#endif /*COLORED_DE_BRUIJN_GRAPH_*/