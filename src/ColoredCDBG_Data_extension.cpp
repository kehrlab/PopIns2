#include "ColoredCDBG_Data_extension.h"



UnitigExtension::UnitigExtension() : ID(0), entropy(-1), dfs_color('w'), dfs_ancestor(0), dfs_discovertime(0), dfs_finishtime(0) {}


void UnitigExtension::join(const UnitigColorMap<UnitigExtension> &um_dest, const UnitigColorMap<UnitigExtension> &um_src){
    // Irrelevant at the moment.
}


void UnitigExtension::sub(UnitigExtension* data_dest, const UnitigColors &uc_dest, const UnitigMapBase &um_dest, const UnitigColorMap<UnitigExtension> &um_src, const bool last_extraction){
    // Irrelevant at the moment.
}


string UnitigExtension::serialize() const{
    string s = std::to_string(getID());
    return s;
}


/*!
 * \fn      uint8_t setNeighborPairFromBases(const uint8_t predecessor, const uint8_t successor)
 * \brief   Putting the bitmasks of predecessor and successor in logical conjunction. The result is a bitmask encoding
 *          both, e.g. 00001011 (predecessor G) & 00001101 (successor C) will be encoded in the 00001001(G+C) pair.
 * \return  an 8-bit int encoding a predecessor-successor pair
 */
uint8_t setNeighborPairFromBases(const uint8_t predecessor, const uint8_t successor){
    return (predecessor & successor);
}


/*!
 * \fn      pair<uint8_t, uint8_t> getBasesFromNeighborPair(const uint8_t neighborPair)
 * \brief   This function separates an encoded predecessor-successor pair. Therefore the pair code has to be put in
 *          logical disjunction with the opposite neighbor's blank bitmask.
 * \return  the predecessor and successor separated (wrapped in a std::pair)
 */
pair<uint8_t, uint8_t> getBasesFromNeighborPair(const uint8_t neighborPair){
    return make_pair<uint8_t, uint8_t>(neighborPair | (0b00000011), neighborPair | (0b00001100));
}
