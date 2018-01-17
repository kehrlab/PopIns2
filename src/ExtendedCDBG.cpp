#include "ExtendedCDBG.h"
#include <unordered_set>



ExtendedCDBG::ExtendedCDBG(int kmer_length, int minimizer_length): CompactedDBG< UnitigExtension >(kmer_length, minimizer_length){
}


void ExtendedCDBG::init_ids(){
    unsigned i=0;
    for (auto& unitig : *this){
        const UnitigExtension ue(i);
        unitig.setData(&ue);
        ++i;
    }
}


void ExtendedCDBG::print_ids(){
    std::cout << "[DEBUG] ";
    for (auto &unitig : *this){
        unsigned id = unitig.getData()->getID();
        std::cout << id << ", ";
    }
    std::cout << std::endl;
}


/*!
 * \fn      size_t ExtendedCDBG::count_connected_components()
 * \brief   This function computes the amount of distict connected components in the graph. The connected_components()
 *          function need to be ran successfully to call count_connected_components(). This function is mainly for
 *          debug and test purposes.
 * \return  number of distinct connected components
 */
size_t ExtendedCDBG::count_connected_components(){
    std:unordered_set<unsigned> unique_set;
    for (auto &unitig : *this){
        unique_set.insert(seqan::findSet(UF, unitig.getData()->getID()));
    }
    return unique_set.size();
}


/*!
 * \fn      bool ExtendedCDBG::connected_components(const CDBG_Build_opt &graph_options)
 * \brief   This function computes the connected components for at the current state of the graph.
 * \ref     seqan/include/seqan/misc/union_find.h
 * \return  0 if successful
 */
bool ExtendedCDBG::connected_components(const CDBG_Build_opt &graph_options){
    // initiate UF structure
    if (graph_options.verbose) std::cout << "[VERBOSE] Initiating UNION-FIND" << std::endl;
    resize(UF, (*this).size());

    // run UF merges
    if (graph_options.verbose) std::cout << "[VERBOSE] Running UNION-FIND" << std::endl;
    for (auto IUnitig = (*this).begin(); IUnitig != (*this).end(); ++IUnitig){
        // TODO: progress indicator here

        /*  Get all successors of a node (unitig).
         *  I could have used ancestors as well but it doesn't matter since every link will
         *  be found if the iterator handles each unitig once.
         */
        ForwardCDBG<UnitigExtension, false> neighbours = IUnitig->getSuccessors();
        for (auto INeighbour = neighbours.begin(); INeighbour != neighbours.end(); ++INeighbour)
            seqan::joinSets(UF, seqan::findSet(UF, IUnitig->getData()->getID()), seqan::findSet(UF, INeighbour->getData()->getID()));
    }

    return 0;
}

