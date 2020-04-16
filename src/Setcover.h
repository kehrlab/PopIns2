/**
 *  @file    src/setcover.h
 *  @brief   This file implements the setcover class
 */
#ifndef SETCOVER_
#define SETCOVER_

#include "UnitigExtension.h"
#include <unordered_set>
#include <unordered_map>

#include "prettyprint.h"        // TODO: delete at release



/**
 *  @class    Setcover
 *  @brief    This class implements a setcover
 */
class Setcover{

    typedef std::unordered_map<unsigned, size_t> path_t;
    typedef std::unordered_set<unsigned> setcover_t;

    size_t _min_kmer_contribution;

    setcover_t _setcover;

    path_t _current_path;

public:
    /**
     *          Constructors
     *  @brief  The Setcover::Setcover determines the smallest possible contig
     *          the PopIns2 merge module can return.
     *  @param  threshold is the minimum amount of novel kmers a path has to
     *          contribute to the setcover to be included.
     */
    Setcover() : _min_kmer_contribution(62) {}      // 62 kmers = 125 bases = 2k-1 (for a default k=63)
    Setcover(size_t threshold) : _min_kmer_contribution(threshold) {}


    /**
     *          Function to add an element to the setcover
     *  @param  id is a unitig id
     *  @param  nb_kmers is the amount of kmers of a unitig
     */
    void add(unsigned id, size_t nb_kmers){_current_path.insert(std::pair<unsigned, size_t>(id, nb_kmers));}


    /**
     *          Function to add an element to the setcover
     *  @param  ucm is a unitig
     */
    void add(const UnitigColorMap<UnitigExtension> &ucm){

        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        this->add(ue->getID(), ucm.size - (ucm.getGraph()->getK() - 1));
    }


    /**
     *          Function to test whether a new path should be included into the setcover
     *  @brief  The function checks whether there is enough novel kmer contribution
     *          of the vurrent path with respect to the setcover. If true then
     *          the elements of the current path are flushed into the setcover.
     *  @return bool; true if path was incorporated into the setcover
     */
    bool test(){

        if(!has_min_contribution()){

            _current_path.clear();

            /* DEBUG */ std::cout << "SETCOVER: Path was NOT included."

            return false;
        }

        for (path_t::const_iterator cit=_current_path.cbegin() ; cit!=_current_path.cend(); ++cit)
            _setcover.insert(cit->first);

        _current_path.clear();

        /* DEBUG */ std::cout << "SETCOVER: Path was included."

        return true;
    }


    /**
     *          This function writes a CSV file of IDs that are included in the setcover
     *  @param  ofile_prefix is the prefix of file to be written
     */
    void write(const std::string ofile_prefix) const{

        std::string fname = ofile_prefix+".setcover.csv";
        ofstream ofile(fname);

        if (ofile.is_open()){

            ofile << "ID,Colour" << endl;

            for (auto cit=_setcover.cbegin(); cit != _setcover.cend(); ++cit)
                ofile << *cit << ",red" << endl;

            ofile.close();
        }
        else cerr << "[popins2 merge] [Setcover::write]: Unable to open output file.";
    }

private:

    void print_current_path() const{prettyprint::print(_current_path);}


    /**
     *          This functions checks if the cumulative sum of kmers that are
     *          not in the setcover yet surpass (GEQ) the minimum threshold.
     *  @return true if path has enough novel contribution
     */
    bool has_min_contribution() const{

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


    /**
     *          This function checks if a given id is in the setcover
     *  @param  id is the id to search for in the setcover
     *  @return true if id is already in setcover
     */
    bool contains(const unsigned id) const{return (_setcover.find(id) != c_.cend()) ? true : false;}

};


#endif /*SETCOVER_*/
