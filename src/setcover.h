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

// NOTE: delete at release
#include "prettyprint.h"


/*!
* \class    Setcover
* \brief    This class implements a setcover for the ccDBG's unitigs.
*/
template <typename TSetcover     = std::unordered_set<unsigned>,
          typename TCurrentNodes = std::vector<unsigned>,
          typename TStartNodes   = std::vector<unsigned> >
class Setcover{
public:
    // Constructors
    Setcover() : t_(2) {}
    Setcover(unsigned t) : t_(t) {}

    /*!
    * \fn       Setcover.hasMinContribution()
    * \brief    Check if the current path contains enough new elements with respect to the set cover
    *           to incorporate the container's elements (current path) into the set cover. Returns true
    *           if the amount of new elements is larger than or equal to the threshold.
    * \return   bool
    */
    bool hasMinContribution() const{
        unsigned novel_elements = 0;
        auto cit = current_path_.cbegin();
        while (cit != current_path_.cend()){
            if (!contains(*cit))
                ++novel_elements;
            if (novel_elements >= t_)
                return true;
            ++cit;
        }
        return false;
    }

    // add/delete/clear functions for the current path container
    void add(unsigned u){current_path_.push_back(u);}
    void add(const UnitigColorMap<UnitigExtension> &ucm){
        DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
        UnitigExtension* ucm_ue = ucm_da->getData(ucm);
        current_path_.push_back(ucm_ue->getID());
    }
    void addStartnode(unsigned s){s_.push_back(s);}
    void del(){current_path_.pop_back();}
    void clear(){current_path_.clear();}

    // Unify the the elements of the current path with the set cover
    void incorporate(){c_.insert(current_path.cbegin(), current_path.cend());}

    void print_CSV(){
        ofstream csv_f("contigs.csv");
        csv_f << "ID,Colour" << endl;
        for (auto it=c_.cbegin(); it != c_.cend(); ++it) csv_f << *it << ",red" << endl;
        for (unsigned n=0; n < s_.size(); n++) csv_f << s_[n] << ",green" << endl;
        csv_f.close();
    }

    void print_current(){prettyprint::print(current_path_);}

private:
    // minimum contribution of a path, checked in Setcover::check()
    unsigned t_;

    // set cover
    TSetcover c_;

    // container for the curent path
    TCurrentNodes current_path_;

    // Container for traversed startnodes
    TStartNodes s_;

    bool contains(const unsigned e) const{
        if (c_.find(e)!=c_.end())
            return true;
        return false;
    }
};


#endif /*SETCOVER_*/