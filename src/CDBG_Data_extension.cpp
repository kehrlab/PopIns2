#include "CDBG_Data_extension.h"



void UnitigExtension::join(const UnitigMap< UnitigExtension >& um_dest, const UnitigMap< UnitigExtension >& um_src){
    if (um_dest.getData()->getConnectedcomponent() < um_src.getData()->getConnectedcomponent())
        um_dest.setData(um_src.getData());
}


void UnitigExtension::sub(const UnitigMap< UnitigExtension >& um_src, UnitigExtension& new_data, const bool last_extraction) const{
    // TODO
}