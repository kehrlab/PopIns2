#include "Setcover.h"



void Setcover::add(const UnitigColorMap<UnitigExtension> &ucm){

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    this->add(ue->getID(), ucm.size - (ucm.getGraph()->getK() - 1));
}


bool Setcover::test(){

    if(!has_min_contribution()){

        _current_path.clear();

        DEBUG_PRINT_STATUS("[popins2 merge][Setcover::test] Path was NOT included.");

        return false;
    }

    for (path_t::const_iterator cit=_current_path.cbegin() ; cit!=_current_path.cend(); ++cit)
        _setcover.insert(cit->first);

    _current_path.clear();

    DEBUG_PRINT_STATUS("[popins2 merge][Setcover::test] Path was included.");

    return true;
}


void Setcover::write(const std::string ofile_prefix) const{

    std::string fname = ofile_prefix+".setcover.csv";
    ofstream ofile(fname);

    if (ofile.is_open()){

        ofile << "ID,Colour" << endl;

        for (auto cit=_setcover.cbegin(); cit != _setcover.cend(); ++cit)
            ofile << *cit << ",red" << endl;

        ofile.close();
    }
    else cerr << "[popins2 merge][Setcover::write]: Unable to open output file.";
}


inline bool Setcover::has_min_contribution() const{

    size_t novel_kmers = 0;

    path_t::const_iterator cit = _current_path.cbegin();
    while (cit != _current_path.cend()){

        if (!contains(cit->first))              // first  = ID
            novel_kmers += cit->second;         // second = #kmers of the unitig with ID

        if (novel_kmers >= _min_kmer_contribution)
            return true;

        ++cit;
    }

    return false;
}


// EOF
