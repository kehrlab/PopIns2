#include "CDBG_Data_extension.h"



UnitigExtension::UnitigExtension(const unsigned initID, const float initEntropy) : ID(initID), entropy(initEntropy) {}


void UnitigExtension::join(const UnitigMap< UnitigExtension >& um_dest, const UnitigMap< UnitigExtension >& um_src){
    /*
     * Irrelevant at the moment.
    */
}


void UnitigExtension::sub(const UnitigMap< UnitigExtension >& um_src, UnitigExtension& new_data, const bool last_extraction) const{
    /*
     * Irrelevant at the moment.
    */
}


string UnitigExtension::serialize() const{
    string s = std::to_string(getID());
    return s;
}