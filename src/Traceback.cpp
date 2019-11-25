#include "Traceback.h"
#include <iostream>
#include <fstream>



void Traceback::join(const Traceback &t){

    for (auto it = t.cbegin(); it != t.cend(); ++it)
        this->push_back(*it);
}


void Traceback::cutconcat(string &s, const VSequences &path, const size_t k) const{
    // since the traceback stored the sequences from sink to source, the concatination has to run backwards over the vector
    for (auto rit = path.crbegin(); rit != path.crend(); ++rit){
        if (rit == path.crbegin())
            s+=(*rit);
        else
            s+=rit->substr(k-1);
    }
}


/*!
 * \fn      bool Traceback::write(const std::string &filename, ofstream &ofs) const
 * \return  true if write was successful
 */
bool Traceback::write(ofstream &ofs, const size_t k, size_t &counter) const{
    if (ofs.is_open()){

#ifdef DEBUG
        prettyprint::print(pathseqs);
#endif // DEBUG

        // loop through all traceback paths found in one node
        for (auto it=cbegin(); it!=cend(); ++it){
            std::string seq;
            cutconcat(seq, *it, k);
            ++counter;

            ofs << ">contig_" << counter << "\n";
            ofs << seq << "\n";
        }
        return 1;
    }
    else
        cout << "Error: Unable to open file." << endl;
    return 0;
}


void Traceback::rearrange(const Traceback &bw, const Traceback &fw){
    //cout << "__________________________" << endl;
    VSequences tmp_v;
    Traceback::const_iterator bw_p = bw.cbegin();   // a valid solution should only contain one vector in VVSequences
    Traceback::const_iterator fw_p = fw.cbegin();   // a valid solution should only contain one vector in VVSequences
    tmp_v.insert(tmp_v.cend(), fw_p->cbegin() , --(fw_p->cend()));
    tmp_v.insert(tmp_v.cend(), bw_p->crbegin(),    bw_p->crend());
    this->push_back(tmp_v);     // tmp_v is actually reverse now but cutconcat() traverses backwards later

    //prettyprint::print(tmp_v);
    //cout << "__________________________" << endl;
}
