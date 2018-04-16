/*!
* \file    src/CompactedDBG_Graph_extension.h
* \brief   Library for a compacted de Bruijn Graph using data extension.
*
*/
#ifndef COMPACTED_DBG_GRAPH_EXTENSION_
#define COMPACTED_DBG_GRAPH_EXTENSION_

#include "CompactedDBG_Data_extension.h"

#include <cstdint>                      /* uint8_t */



// =========================
// Global functions
// =========================


// =========================
// Struct
// =========================
/*!
* \class        ExtendedCDBG
* \headerfile   src/CompactedDBG_Graph_extension.h
* \brief        Struct to store a compacted DBG plus data extension.
*/
struct ExtendedCDBG : public CompactedDBG<DataExtension> {

    public:
        // constructor
        ExtendedCDBG(int kmer_length, int minimizer_length) :   CompactedDBG<DataExtension,void> (kmer_length, minimizer_length),
                                                                BFS_MAX_DIST(kmer_length<<1),
                                                                id_init_status(false)
                                                                {}

        void init_ids();

        uint8_t whereToGo(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src) const;

        void clear_traversal_attributes();

        bool small_bubble_removal();
        bool BFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction);
        bool BFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist);

    private:

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        const size_t BFS_MAX_DIST;   // default: two times k-mer length

        bool id_init_status;

};





















#endif /*COMPACTED_DBG_GRAPH_EXTENSION_*/