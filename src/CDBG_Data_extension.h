/*!
* \file    src/CDBG_Data_extension.h
* \brief   Library for extending the functionality of the CompactedDBG by implementing the data wrapper CDBG_Data_t
*
*/
#ifndef UNITIG_EXTENSION_H_
#define UNITIG_EXTENSION_H_


// =========================
// Includes + defines
// =========================
#include <bifrost/CompactedDBG.hpp>     /* has the CDBG_Data_t template */
#include <string>                       /* std::to_string() */
#include <tuple>                        /* std::Pair */



// =========================
// Structs
// =========================
/*!
* \class        UnitigExtension
* \headerfile   src/CDBG_Data_extension.h
* \brief        Struct extending the functionality of the CompactedDBG.
* \details      The struct has to inherit from the CDBG_Data_t struct and implement at least all its virtual functions.
* \ref          https://github.com/pmelsted/bfgraph/blob/master/src/CompactedDBG.hpp
*/
struct UnitigExtension : public CDBG_Data_t<UnitigExtension> {

    private:
        unsigned ID;

        float entropy;

        /* A neighbourPair (uint8_t) is a 8-bit int storing the predecessor in
         * the third and fourth last bits (xxxxXX11) and the successor in the last
         * and second last bits (xxxx11XX), where the two bits encode A, C, G, and T.
         * The bits of the partner neighbor is always set to 11, so it is consistent
         * with the bit operations defined for neighbor pairs.
         */
        std::vector<uint8_t> neighborPairs;

    public:
        UnitigExtension(const unsigned initID = 0, const float initEntropy = -1);

        unsigned getID() const {return ID;}
        void setID(const unsigned id) {ID = id;}

        float getEntropy() const {return entropy;}
        void setEntropy(const float e) {entropy = e;}

        void join(const UnitigMap<UnitigExtension>& um_dest, const UnitigMap<UnitigExtension>& um_src);
        void sub(const UnitigMap<UnitigExtension>& um_src, UnitigExtension& new_data, const bool last_extraction) const;
        string serialize() const;

};


// =========================
// Global functions
// =========================
uint8_t setNeighborPairFromBases(const uint8_t predecessor, const uint8_t successor);
std::pair<uint8_t, uint8_t> getBasesFromNeighborPair(const uint8_t neighborPair);



#endif /*UNITIG_EXTENSION_H_*/