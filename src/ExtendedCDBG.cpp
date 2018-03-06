#include "ExtendedCDBG.h"
#include <unordered_set>
#include <unordered_map>



// default constructor
ExtendedCDBG::ExtendedCDBG(int kmer_length, int minimizer_length): CompactedDBG< UnitigExtension >(kmer_length, minimizer_length){
    init_status = false;    // IDs are not initiated at construction time (see init_ids())
}


void ExtendedCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto& unitig : *this){
        UnitigExtension* ue = unitig.getData();      // ue is a POINTER to a UnitigExtension
        ue->setID(i);
        ++i;
    }
    init_status = true;
}


void ExtendedCDBG::print_ids(){
    if (is_init() == true){
        std::cout << "[DEBUG] ";
        for (auto &unitig : *this){
            unsigned id = unitig.getData()->getID();
            std::cout << id << ", ";
        }
        std::cout << std::endl;
    }
    else
        cerr << "[WARNING] Unitig IDs were not printed because they are not initialized." << endl;
}


bool ExtendedCDBG::is_init(){
    return init_status;
}


/*!
 * \fn      size_t ExtendedCDBG::count_connected_components()
 * \brief   This function computes the amount of distict connected components in the graph. The connected_components()
 *          function need to be ran successfully to call count_connected_components(). This function is mainly for
 *          debug and test purposes.
 * \return  number of distinct connected components
 */
size_t ExtendedCDBG::count_connected_components(){
    std::unordered_set<unsigned> unique_set;
    for (auto &unitig : *this){
        unique_set.insert(seqan::findSet(UF, unitig.getData()->getID()));
    }
    return unique_set.size();
}


/*!
 * \fn      bool ExtendedCDBG::connected_components(const CDBG_Build_opt &graph_options)
 * \brief   This function computes the connected components for at the current state of the graph.
 * \ref     seqan/include/seqan/misc/union_find.h
 * \return  true if successful
 */
bool ExtendedCDBG::connected_components(const CDBG_Build_opt &graph_options){
    // initiate UF structure
    if (graph_options.verbose) std::cout << "[VERBOSE] Initiating UNION-FIND" << std::endl;
    resize(UF, (*this).size());
    //std::cout << "UF size " << seqan::length(UF._values) << std::endl;    //[DEBUG]

    // run UF merges
    if (graph_options.verbose) std::cout << "[VERBOSE] Running UNION-FIND" << std::endl;
    for (auto &unitig : *this){
        // TODO: progress indicator here

        /*  Get all predecessors AND successors of a node (unitig).
         *  I need to iterate through both since either one of them could miss
         *  links, and therefore split components, where both unitigs
         *  "face each other", e.g.:
         *        u1 ----------->
         *                 <-------------- u2
         * which is in GFA:
         *      L   u1  +   u2  -
         *      L   u2  +   u1  -
         */
        size_t unitig_id = unitig.getData()->getID();

        BackwardCDBG<UnitigExtension, false> predecessors = unitig.getPredecessors();
        for (auto &it_pre : predecessors){
            size_t pre_id = it_pre.getData()->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, pre_id));
        }

        ForwardCDBG<UnitigExtension, false> successorss = unitig.getSuccessors();
        for (auto &it_suc : successorss){
            size_t suc_id = it_suc.getData()->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, suc_id));
        }
    }

    return true;
}


void ExtendedCDBG::print_components(){
    for (auto &unitig : *this){
        size_t uid = unitig.getData()->getID();
        unsigned ucc = seqan::findSet(UF, uid);
        cout << "ID:" << uid << " | " << "CC:" << ucc << endl;
    }
}



/*!
 * \fn          float ExtendedCDBG::entropy(const std::string &sequence)
 * \brief       This function computes an entropy for a given string that can be used to filter/mark low complexity
 *              sequences. If all dimers are equaly distributed the entropy is high ("highly chaotic system"), if
 *              all dimers follow a certain pattern the entropy is low ("highly ordered system"). We'd probably like
 *              to mark low entropy unitigs since they have a chance to disrupt/branch the de Bruijn Graph.
 * \remark      Function taken from PopIns.
 * \return      The entropy [0,1] of all binucleotides.
 */
float ExtendedCDBG::entropy(const std::string &sequence){
    // create a dictionary counting the occurrence of all dinucleotides
    unordered_map<std::string, unsigned> diCounts(16);
    unsigned counted = 0;
    for (unsigned i = 0; i < sequence.length()-1; ++i){
        std::string dimer = sequence.substr(i,2);
        if (sequence[i]!='N' && sequence[i+1]!='N'){
            // set if dimer not in counter table yet
            if(diCounts.find(dimer) == diCounts.end()){
                diCounts[dimer] = 1;
                counted++;
            }
            // otherwise increase
            else{
                diCounts[dimer] += 1;
                counted++;
            }
        }
    }

    // calculate the entropy for dinucleotide counts
    unsigned entropy = 0;
    for(unordered_map<std::string,unsigned>::const_iterator it = diCounts.cbegin(); it != diCounts.cend(); ++it){
        if (it->second == 0) continue;
        float p = float(it->second) / counted;
        entropy -= p * log(p) / log(2);
    }

    return entropy / 4;
}






