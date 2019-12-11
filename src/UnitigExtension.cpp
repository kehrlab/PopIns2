#include "UnitigExtension.h"


/*
inline void UnitigExtension::clear(const UnitigColorMap<UnitigExtension> &um_dest){
    // defaults are fine
}


inline void UnitigExtension::concat(const UnitigColorMap< UnitigExtension >& um_dest,
                                    const UnitigColorMap< UnitigExtension >& um_src){
    // defaults are fine
}


inline void UnitigExtension::merge(const UnitigColorMap< UnitigExtension >& um_dest,
                                   const const_UnitigColorMap< UnitigExtension >& um_src){
    // defaults are fine
}


inline void UnitigExtension::extract(const UnitigColors* uc_dest,
                                     const UnitigColorMap< UnitigExtension >& um_src,
                                     const bool last_extraction){
    // defaults are fine
}
*/


/**/
string UnitigExtension::serialize(const const_UnitigColorMap< UnitigExtension >& um_src) const{

    const DataAccessor<UnitigExtension>* da = um_src.getData();
    const UnitigExtension* ue = da->getData(um_src);

    string field = "ID:i:" + std::to_string(ue->getID()) + "\tEN:f:" + std::to_string(ue->getEntropy());

    return field;
}
/**/


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
    return make_pair<uint8_t, uint8_t>(neighborPair | (0x03), neighborPair | (0x0C));
}
