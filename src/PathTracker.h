/*!
* \file     src/PathTracker.h
* \brief    This file contains the class PathTracker
*/
#ifndef PATHTRACKER_
#define PATHTRACKER_

#include "UnitigExtension.h"
#include "prettyprint.h"        // NOTE: delete at release



class PathTracker{
public:
    //============
    // FUNCTIONS
    //============
    PathTracker(size_t s) : nb_samples(s) {kmers_per_color_.resize(s, 0);}

    void clear(){
        nb_path_kmers_ = 0;
        nb_unitigs_ = 0;
        avg_entropy_ = 0;
        reset_kmer_counter_vector();
    }

    /*!
     * \fn      PathTracker.insert()
     * \brief   incorporates the color and sequence information of a new UnitigColorMap into the PathTracker instance
     */
    void insert(const const_UnitigColorMap<UnitigExtension> &ucm);

    float getEntropy() const {return avg_entropy_;}

private:
    //============
    // MEMBER
    //============
    // this vector stores for all colors the relative amount of kmers in the current path
    std::vector<unsigned> kmers_per_color_;

    unsigned nb_path_kmers_ = 0;

    size_t nb_samples;

    float avg_entropy_ = 0;

    unsigned nb_unitigs_ = 0;

    //============
    // FUNCTIONS
    //============
    void reset_kmer_counter_vector(){for (auto &i : kmers_per_color_) i=0; }

    /*!
     * \fn      PathTracker.maxColorRatio()
     * \brief   highest color ratio (most present color in path)
     */
    float maxColorRatio() const{
        float maxRatio = 0;
        for (auto i : kmers_per_color_)
            if (i/(float)nb_path_kmers_ > maxRatio)
                maxRatio = i/(float)nb_path_kmers_;
        return maxRatio;
    }

    /*!
     * \fn      PathTracker.minColorRatio()
     * \brief   lowest non-zero color ratio (lowest present color in path)
     */
    float minColorRatio() const{
        float minRatio = 1;
        for (auto i : kmers_per_color_)
            if (i/(float)nb_path_kmers_ < minRatio && i/(float)nb_path_kmers_ > 0)
                minRatio = i/(float)nb_path_kmers_;
        return minRatio;
    }

    unsigned nb_colorsInPath() const {unsigned counter=0; for (auto &i : kmers_per_color_) if (i!=0) ++counter; return counter;}
};





#endif /*PATHTRACKER_*/