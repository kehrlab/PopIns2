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
template <typename TContainer = std::set>
class Setcover{
    using iterator = TContainer::iterator;
    using const_iterator = TContainer::const_iterator;
    using key_type = TContainer::key_type;

public:
    Setcover(unsigned p) : p_(p) {}

private:
    unsigned p_ = 2;
    TContainer c_;

    void insert(key_type e){c_.insert(e);}
    template <typename TCopyContainer> void insert_range(TCopyContainer::iterator &begin, TCopyContainer::iterator &end){c_.insert(begin, end);}
    template <typename TCopyContainer> void insert_range(TCopyContainer::const_iterator &cbegin, TCopyContainer::const_iterator &cend){c_.insert(cbegin, cend);}
    
    template <typename TCheckContainer> bool check();

};


#endif /*SETCOVER_*/