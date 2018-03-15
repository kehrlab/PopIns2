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


// =========================
// Structs
// =========================
/*!
* \class        ExtendedCDBG
* \headerfile   src/ExtendedCDBG.h
* \brief        Struct to store a CompactedDBG plus data extensions.
*/
struct ExtendedCDBG : public CompactedDBG<UnitigExtension> {

    private:

        bool init_status;
        bool isKmerCovInit;

        seqan::UnionFind<unsigned> UF;

        size_t dfs_time;
        bool dfs_passed;


    public:
        ExtendedCDBG(int kmer_length = 31, int minimizer_length = 23);      // hidden inits! (see definition)

        void print_unitig_info();

        void init_ids();
        void print_ids();
        bool is_init() const {return init_status;}

        bool connected_components(const CDBG_Build_opt &graph_options);
        size_t count_connected_components();
        seqan::UnionFind<unsigned> getUF() const {return UF;}

        float entropy(const std::string &sequence);

        void dfs(UnitigMap<UnitigExtension> &um);
        void dfs_visit(UnitigMap<UnitigExtension> &um);
        bool is_dfs_passed() const {return dfs_passed;}

        void traceback(vector<unsigned> &vec, UnitigMap<UnitigExtension> &um_sink);

        bool annotate_kmer_coverage(const vector<string> &sample_fastx_names);

        void init_kmer_cov();
};



#endif /*EXTENDED_CDBG_H_*/