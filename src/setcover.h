/*!
* \file    src/setcover.h
* \brief   This file separates the setcover object from the main code
*
*/
#ifndef SETCOVER_
#define SETCOVER_

#include <set>


/*!
* \class    Setcover
* \brief    This class is a wrapper for a container to sto
*/
template <typename TContainer = std::set<unsigned> >
class Setcover{
    using iterator = TContainer::iterator;
    using const_iterator = TContainer::const_iterator;
    using key_type = TContainer::key_type;

public:
    Setcover(unsigned t) : t_(t) {}

    /*!
    * \fn       Setcover.check()
    * \brief    Check if an input container contains enough new elements, with respect to the set cover,
    *           to unify the container with the set cover.
    * \return   bool; true if #new elements is >= threshold
    */
    template <typename TCheckContainer> bool check(TCheckContainer::const_iterator &cbegin, TCheckContainer::const_iterator &cend){
        unsigned novel_elements = 0;
        while (cbegin != cend){
            if (c_.find(*cbegin) != c_.end())
                ++novel_elements;
            if (novel_elements >= this->t)
                return true;
            ++cbegin;
        }
        return false;
    }

private:
    unsigned t_ = 2;
    TContainer c_;

    /* Insert an element or a range of elements into the set cover. */
    void insert(key_type e){c_.insert(e);}
    template <typename TInsertContainer> void insert_range(TInsertContainer::const_iterator &cbegin, TInsertContainer::const_iterator &cend){c_.insert(cbegin, cend);}

};


#endif /*SETCOVER_*/