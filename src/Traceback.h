/*!
* @file    src/Traceback.h
* @brief   Container class for managing the traceback of the dBG traversal
*
*/
#ifndef TRACEBACK_
#define TRACEBACK_

#include <iostream>
#include <fstream>

using namespace std;



// =========================
// Structs
// =========================

/*!
* @class        Traceback
* @headerfile   src/Traceback.h
* @brief        Class to manage the metadata for the DFS traceback.
*/
class Traceback{

    typedef uint8_t direction_t;

private:
    // ----------
    // | member |
    // ----------

    string _contig;                         // string vectors for paths

    const direction_t _d;                   // direction of the traversal

    const size_t _k;                        // kmer length of the dBG


public:
    // ---------------
    // | constructor |
    // ---------------

    Traceback(const direction_t d, const size_t k) : _d(d), _k(k) {}

    // -------------
    // | functions |
    // -------------

    /**
     *  Function to trim and concatenate the current unitig with the final contig
     *  @brief  This function concatenates the path from sink to source.
     *  @param  unitig is the string to be concatenated with the final contig
     *  @param  startnode is an indicator whether unitig was the first node of the traversal
     *          and therefore should not be trimmed.
    */
    void add(const string &unitig, const bool startnode = false);

    /**
    *   Function to add 'N's to the final contig
    *   @brief  This function adds 'N' bases to the contig to mask an unknown sequence fragment.
    */
    void addN();

    /**
    *   Function to add the full sequence of a given string
    *   @brief  This function adds the full sequence string without trimming the kmer overlap.
    *           This function covers the edge case that a sink (unitig) of a path can directly
    *           follow a jump through a LECC. In this case trimming k-1 bases off the sink would
    *           destroy valuable (for the linking module might crucial) information.
    */
    void addFullSink(const string &unitig);

    /**
     *  Function to write the current contig
     *  @param  of is a file stream to write the contig into
    */
    void write(ofstream &of, const size_t counter) const;

    string print(){cout << _contig << endl; return _contig;}    // DEBUG
};










#endif /*TRACEBACK_*/
