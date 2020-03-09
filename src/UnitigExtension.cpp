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
