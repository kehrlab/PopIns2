/*!
* \file    src/setcover.h
* \brief   This file separates the setcover object from the main code
*
*/
#ifndef SETCOVER_
#define SETCOVER_

#include "UnitigExtension.h"
#include <unordered_set>
#include <vector>

#include <seqan/basic.h>

// NOTE: delete at release
#include "prettyprint.h"



/*!
* \class    Setcover
* \brief    This class implements a setcover for the ccDBG's unitigs.
*/
template <typename TSetcover     = std::unordered_set<unsigned>,
          typename TCurrentPath  = std::unordered_map<unsigned, unsigned> >
class Setcover{
public:
    // Constructors
    Setcover() : t_(255) {}
    Setcover(unsigned t) : t_(t) {}

    /**
     *          Checks if the cumulative sum of unseen kmers of the current path p is LEQ a threshold.
     * @return  true if cumsum(p) is larger than or equal threshold t_
     */
    bool hasMinContribution(bool verbose = false) const{
        unsigned novel_kmers = 0;
        auto cit = current_path_.cbegin();
        while (cit != current_path_.cend()){
            if (!contains(cit->first)){
                novel_kmers += cit->second;
                if (verbose) cout << cit->first << " is novel" << endl;
            }
            else{
                if (verbose) cout << cit->first << " is not novel" << endl;
            }
            if (novel_kmers >= t_){
                if (verbose) cout << "Found enough novel kmers" << endl;
                return true;
            }
            ++cit;
        }
        if (verbose) cout << "Found not enough novel kmers" << endl;
        return false;
    }

    // add/delete/clear functions for the current path container
    void add(unsigned err){(void)err; SEQAN_THROW(seqan::RuntimeError("Argument error in sc.add() "));}
    void add(unsigned id, unsigned nb_kmers){current_path_.insert(std::pair<unsigned, unsigned>(id, nb_kmers));}
    void add(const UnitigColorMap<UnitigExtension> &ucm){
        DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
        UnitigExtension* ucm_ue = ucm_da->getData(ucm);
        add(ucm_ue->getID(), ucm.len);
    }
    void del(unsigned id){current_path_.erase(id);}
    void del(const UnitigColorMap<UnitigExtension> &ucm){
        DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
        UnitigExtension* ucm_ue = ucm_da->getData(ucm);
        del(ucm_ue->getID());
    }
    void clear(){current_path_.clear();}

    bool empty(){return c_.empty();}

    // Unify the the elements of the current path with the set cover
    void incorporate(bool verbose = false){
        typename TCurrentPath::const_iterator cit = current_path_.cbegin();
        for ( ; cit != current_path_.cend(); ++cit) c_.insert(cit->first);
        if (verbose) cout << "Inserted kmers of following novel elements " << endl;
        if (verbose) print_current();
        if (verbose) cout << endl;
    }

    void write_CSV(const std::string s) const{
        std::string csv_file_name = s+".csv";
        ofstream csv_f(csv_file_name);
        csv_f << "ID,Colour" << endl;
        for (auto it=c_.cbegin(); it != c_.cend(); ++it) csv_f << *it << ",red" << endl;
        csv_f.close();
    }

    void print_current(){prettyprint::print(current_path_);}

    bool contains(const unsigned id) const{
        if (c_.find(id)!=c_.end())
            return true;
        return false;
    }

private:
    // minimum contribution of a path, checked in Setcover::check()
    unsigned t_;

    // set cover
    TSetcover c_;

    // container for the curent path
    TCurrentPath current_path_;

};



/*!
* \class    Setcover_OLD
* \brief    This class implements a setcover for the ccDBG's unitigs.
*/
template <typename TSetcover     = std::unordered_set<unsigned>,
          typename TCurrentPath  = std::unordered_set<unsigned> >
class Setcover_OLD{
public:
    // Constructors
    Setcover_OLD() : t_(2) {}
    Setcover_OLD(unsigned t) : t_(t) {}

    /**
     *          Checks if the current path p contains enough new elements with respect to the set cover c
     *          to incorporate the elements of p into the c.
     * @return  true if number of new elements is larger than or equal threshold t_
     */
    bool hasMinContribution(bool verbose = true) const{
        unsigned novel_elements = 0;
        auto cit = current_path_.cbegin();
        while (cit != current_path_.cend()){
            if (!contains(*cit)){
                ++novel_elements;
                if (verbose) cout << *cit << " is novel" << endl;
            }
            else{
                if (verbose) cout << *cit << " is not novel" << endl;
            }
            if (novel_elements >= t_){
                if (verbose) cout << "Found enough novel elements" << endl;
                return true;
            }
            ++cit;
        }
        if (verbose) cout << novel_elements << " of " << current_path_.size() << " are novel to the SC" << endl;
        return false;
    }

    // add/delete/clear functions for the current path container
    void add(unsigned u){current_path_.insert(u);}
    void add(const UnitigColorMap<UnitigExtension> &ucm){
        DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
        UnitigExtension* ucm_ue = ucm_da->getData(ucm);
        current_path_.insert(ucm_ue->getID());
    }
    void del(unsigned u){current_path_.erase(u);}
    void clear(){current_path_.clear();}

    bool empty(){return c_.empty();}

    // Unify the the elements of the current path with the set cover
    void incorporate(bool verbose = true){
        c_.insert(current_path_.cbegin(), current_path_.cend());
        if (verbose) cout << "Inserted following novel elements ";
        if (verbose) print_current();
        if (verbose) cout << endl;
    }

    void write_CSV(){
        ofstream csv_f("contigs.csv");
        csv_f << "ID,Colour" << endl;
        for (auto it=c_.cbegin(); it != c_.cend(); ++it) csv_f << *it << ",red" << endl;
        csv_f.close();
    }

    void print_current(){prettyprint::print(current_path_);}

    bool contains(const unsigned e) const{
        if (c_.find(e)!=c_.end())
            return true;
        return false;
    }

private:
    // minimum contribution of a path, checked in Setcover::check()
    unsigned t_;

    // set cover
    TSetcover c_;

    // container for the curent path
    TCurrentPath current_path_;

};


#endif /*SETCOVER_*/
