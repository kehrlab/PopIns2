/*!
* \file     src/ColorTracker.h
* \brief    This file contains the class ColorTracker
*/
#ifndef COLORTRACKER_
#define COLORTRACKER_

#include "UnitigExtension.h"
#include "prettyprint.h"        // NOTE: delete at release



class ColorTracker{
public:
    ColorTracker(size_t s) : nb_samples(s) {kmer_counts_.resize(s, 0);}

    void clear(){
        total_kmers_ = 0;
        reset();
    }

    /*!
     * \fn      ColorTracker.update()
     * \brief   update() incorporates the kmer colors of a new UnitigColorMap into the reference color counts
     */
    void update(const const_UnitigColorMap<UnitigExtension> &ucm){
        total_kmers_ += ucm.len;

        const UnitigColors* unitig_colors = ucm.getData()->getUnitigColors(ucm);
        size_t max_color_index = unitig_colors->colorMax(ucm);

        // sum up #kmers per color
        for (auto color_id=0; color_id<=max_color_index; ++color_id){
            size_t nb_kmers = unitig_colors->size(ucm, color_id);
            kmer_counts_[color_id] += nb_kmers;
        }
    }

    float maxRatio() const{
        float maxRatio = 0;
        for (auto i : kmer_counts_)
            if (i/(float)total_kmers_ > maxRatio)
                maxRatio = i/(float)total_kmers_;
        return maxRatio;
    }

    float minRatio() const{
        float minRatio = 1;
        for (auto i : kmer_counts_)
            if (i/(float)total_kmers_ < minRatio && i/(float)total_kmers_ > 0)
                minRatio = i/(float)total_kmers_;
        return minRatio;
    }

private:
    // this vector stores for all colors the relative amount of kmers in the current path
    std::vector<unsigned> kmer_counts_;

    // this is the sum of kmers of all unitigs in the current path
    unsigned total_kmers_ = 0;

    size_t nb_samples;

    void reset(){for (auto i : kmer_counts_) i=0; }
};





#endif /*COLORTRACKER_*/