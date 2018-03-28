/*!
* \file    src/ExtendedCDBG.h
* \brief   Library for a compacted de Bruijn Graph using data extension.
*
*/
#ifndef EXTENDED_CDBG_H_
#define EXTENDED_CDBG_H_


// =========================
// Includes
// =========================
#include <seqan/seq_io.h>
#include <seqan/misc/union_find.h>

#include "CDBG_Data_extension.h"

#include "argument_parsing.h"

typedef std::vector<std::vector<UnitigColorMap<UnitigExtension> > > PathSet;
typedef std::vector<UnitigColorMap<UnitigExtension> > UnitigPath;



// =========================
// Structs
// =========================
/*!
* \class        ExtendedCDBG
* \headerfile   src/ExtendedCDBG.h
* \brief        Struct to store a CompactedDBG plus data extensions.
*/
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    private:

        bool init_status;
        bool isKmerCovInit;

        seqan::UnionFind<unsigned> UF;

        size_t dfs_time;
        bool dfs_passed;


    public:
        ExtendedCCDBG(int kmer_length = 31, int minimizer_length = 23);      // hidden inits! (see definition)

        // wrapper to skip the dereferencing via the DataAccessor
        //const UnitigExtension* getDataDirectly() const;
        //UnitigExtension* getDataDirectly(UnitigColorMap &ucm);      // inline

        void print_unitig_info();

        void init_ids();
        void print_ids();
        bool is_init() const {return init_status;}

        bool connected_components(const CCDBG_Build_opt &graph_options);
        size_t count_connected_components();
        seqan::UnionFind<unsigned> getUF() const {return UF;}

        float entropy(const std::string &sequence);

        void dfs(const UnitigColorMap<UnitigExtension> &um);
        void dfs_visit(const UnitigColorMap<UnitigExtension> &um);
        bool is_dfs_passed() const {return dfs_passed;}
        void dfs_traceback(vector<unsigned> &vec, const UnitigColorMap<UnitigExtension> &um_sink);

        void init_kmer_cov();
        bool annotate_kmer_coverage(const vector<string> &sample_fastx_names);

        void clear_path_search_attributes();

        void small_bubble_removal();    // bubble popping main
        bool bfs_with_max_dist(const UnitigColorMap<UnitigExtension> &um, PathSet &pathset, const size_t max_dist);       // inline function to detect bubbles smaller than delta_k
        bool get_reverse_bfs_paths(const UnitigColorMap<UnitigExtension>& um, UnitigPath &up);   // get all paths from sink to source
};








#endif /*EXTENDED_CDBG_H_*/