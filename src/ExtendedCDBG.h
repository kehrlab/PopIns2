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

        seqan::UnionFind<unsigned> UF;


    public:
        ExtendedCDBG(int kmer_length = 31, int minimizer_length = 23);

        void init_ids();
        void print_ids();
        bool is_init();

        bool connected_components(const CDBG_Build_opt &graph_options);
        size_t count_connected_components();
        void print_components();
        seqan::UnionFind<unsigned> getUF() const {return UF;}

        float entropy(const std::string &sequence);
};



#endif /*EXTENDED_CDBG_H_*/