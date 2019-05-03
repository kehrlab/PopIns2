/*!
* \file    src/setcover.h
* \brief   This file separates the setcover object from the main code
*
*/
#ifndef SETCOVER_
#define SETCOVER_

#include "UnitigExtension.h"
#include <set>
#include <vector>


/*!
* \class    Setcover
* \brief    This class implements a setcover for the ccDBG's unitigs.
* \param    
*/
template <typename TSetcoverContainer = std::set<unsigned>, 
          typename TPathContainer     = std::vector<unsigned> >
class Setcover{
public:
    // Constructors
    Setcover() : t_(2) {}
    Setcover(unsigned t) : t_(t) {}

    /*!
    * \fn       Setcover.check()
    * \brief    Check if the current path contains enough new elements, with respect to the set cover,
    *           to unify the container (of current path) with the set cover.
    * \return   bool; true if #new elements is >= threshold
    */
    bool check() const{
        unsigned novel_elements = 0;
        auto cit = current_path.cbegin();

        while (cit != current_path.cend()){
            if (c_.find(*cit) != c_.end())
                ++novel_elements;
            if (novel_elements >= this->t_)
                return true;
            ++cit;
        }
        return false;
    }

    // add/delete/clear functions for the current path container
    void add(unsigned u){current_path.push_back(u);}
    void add(const UnitigColorMap<UnitigExtension> &ucm){
        DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
        UnitigExtension* ucm_ue = ucm_da->getData(ucm);
        current_path.push_back(ucm_ue->getID());
    }
    void del(){current_path.pop_back();}
    void clear(){current_path.clear();}

    // Unify the the elements of the current path with the set cover
    void unify(){c_.insert(current_path.cbegin(), current_path.cend());}

private:
    // minimum contribution of a path, checked in Setcover::check()
    unsigned t_;

    // set cover
    TSetcoverContainer c_;

    // container for the curent path
    TPathContainer current_path;

};


#endif /*SETCOVER_*/