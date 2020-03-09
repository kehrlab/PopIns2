#include "PathTracker.h"



void PathTracker::insert(const const_UnitigColorMap< UnitigExtension >& ucm){
    nb_path_kmers_ += ucm.len;
    ++nb_unitigs_;

    // insert #kmers per color
    const UnitigColors* unitig_colors = ucm.getData()->getUnitigColors(ucm);
    size_t max_color_index = unitig_colors->colorMax(ucm);

    for (size_t color_id=0; color_id<=max_color_index; ++color_id){
        size_t nb_kmers = unitig_colors->size(ucm, color_id);
        kmers_per_color_[color_id] += nb_kmers;
    }

    // update entropy of the current path
    const DataAccessor<UnitigExtension>* da = ucm.getData();
    const UnitigExtension* ue = da->getData(ucm);

    avg_entropy_ = avg_entropy_*((nb_path_kmers_-ucm.len)/nb_path_kmers_) + ue->getEntropy()*(ucm.len/nb_path_kmers_);  // weighted mean
}