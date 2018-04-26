#include "CompactedDBG_Data_extension.h"



// Join method for CompactedDBG
void DataExtension::join(const UnitigMap<DataExtension>& um_dest, const UnitigMap<DataExtension>& um_src){

    // When joining the unitig matching um_src to the unitig matching um_dest,
    // we set um_dest to "not visited" because it will be a new unitig in the graph.

    DataExtension* de_dest = um_dest.getData();

    de_dest->set_not_seen_visited(); // Set the unitig to "not visited"

    DataExtension* de_src = um_src.getData();

    cout << "DEST ID " << de_dest->getID() << " | SRC ID " << de_src->getID() << endl;

    de_dest->setID(de_src->getID());

    cout << "DEST ID " << de_dest->getID() << endl;
}


// Sub method for CompactedDBG
void DataExtension::sub(DataExtension* data_dest, const UnitigMap<DataExtension>& um_src, bool last_extraction){

    // This function creates a new unitig which is a sub-unitig from um_src
    // The new unitig created is set to "not visited" as a measure of precaution
    // (it is already initialed by default to "not visited" in the constructor)

    data_dest->set_not_seen_visited();

    // WARNING: handle is missing for IDs! Not sure yet if needed
    cout << "IS THIS EVER EXECUTED!?" << endl;
}


string DataExtension::serialize() const{
    string s = std::to_string(getID());
    return s;
}

