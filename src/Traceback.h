/*!
* \file    src/Traceback.h
* \brief   Container class for managing the traceback of the dBG traversal
*
*/
#ifndef TRACEBACK_
#define TRACEBACK_

#include <vector>
#include "prettyprint.h"    // TODO: delete at release

using namespace std;

typedef std::vector<std::vector<std::string> >VVSequences;
typedef std::vector<std::string> VSequences;


// =========================
// Structs
// =========================

/*!
* \class        Traceback
* \headerfile   src/ColoredDeBruijnGraph.h
* \brief        Struct to manage the metadata for the DFS traceback.
*/
class Traceback{

private:
    void cutconcat(string &s, const VSequences &path, const size_t k) const;

    VVSequences pathseqs;

public:
    using iterator               = VVSequences::iterator;
    using const_iterator         = VVSequences::const_iterator;
    using const_reverse_iterator = VVSequences::const_reverse_iterator;

    bool recursive_return_status = false;

    bool write(ofstream &ofs, const size_t k, size_t &counter) const;

    void join(const Traceback &t);

    void rearrange(const Traceback &bw, const Traceback &fw);

    bool empty() const {return pathseqs.empty();}
    void push_back(const VSequences &ps) {pathseqs.push_back(ps);}
    const_iterator cbegin() const { return pathseqs.cbegin(); }
    const_iterator cend() const { return pathseqs.cend(); }
    iterator begin() { return pathseqs.begin(); }
    iterator end() { return pathseqs.end(); }
    const_reverse_iterator crbegin() const {return pathseqs.crbegin();}
    const_reverse_iterator crend() const {return pathseqs.crend();}
};










#endif /*TRACEBACK_*/