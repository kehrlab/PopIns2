#include "LECC_Finder.h"



unsigned LECC_Finder::annotate(){
    // sanity check
    if(!this->isGraphInit()){
        std::cerr << "[ERROR][LECC_Finder]: Graph must be initialized with unitig IDs and entropies." << '\n';
        return EXIT_FAILURE;
    }

    unsigned LECC_ = 0;

    for (auto &ucm : *g_){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        // skip high entropy unitigs
        if(ue->getEntropy() >= this->threshold_)
            continue;

        // skip if already LECC-annotated (0 is default LECC identifier)
        if(ue->getLECC() != 0)
            continue;

        // found new LECC
        ++LECC_;
        ue->setLECC(LECC_);

        // traverse predecessors
        for (auto pre : ucm.getPredecessors()){
            this->annotate_recursion(pre, LECC_);
        }

        // traverse successors
        for (auto suc : ucm.getSuccessors()){
            this->annotate_recursion(suc, LECC_);
        }

    }

    return LECC_;
}


void LECC_Finder::annotate_recursion(UnitigColorMap<UnitigExtension> &ucm, const unsigned LECC__){
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    // skip high entropy unitigs
    if(ue->getEntropy() >= this->threshold_)
        return;

    // skip if already LECC-annotated (0 is default LECC identifier)
    if(ue->getLECC() != 0)
        return;

    // found new LECC
    ue->setLECC(LECC__);

    // traverse predecessors
    for (auto pre : ucm.getPredecessors()){
        this->annotate_recursion(pre, LECC__);
    }

    // traverse successors
    for (auto suc : ucm.getSuccessors()){
        this->annotate_recursion(suc, LECC__);
    }
}
