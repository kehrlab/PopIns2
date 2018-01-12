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
#include "CDBG_Data_extension.h"


// =========================
// Structs
// =========================
/*!
* \class        ExtendedCDBG
* \headerfile   src/ExtendedCDBG.h
* \brief        Struct to store a CompactedDBG plus data extensions.
*/
struct ExtendedCDBG : public CompactedDBG<UnitigExtension> {

    public:
        ExtendedCDBG(int kmer_length = 31, int minimizer_length = 23);



};



#endif /*EXTENDED_CDBG_H_*/