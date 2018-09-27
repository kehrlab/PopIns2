#include "ColoredDeBruijnGraph.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>




// default constructor
ExtendedCCDBG::ExtendedCCDBG(int kmer_length, int minimizer_length) :   ColoredCDBG<UnitigExtension> (kmer_length, minimizer_length),
                                                                        id_init_status(false) {
    /* 1) IDs are not initiated at construction time (see init_ids())
     * 2) The UnionFind vector will be empty an construction time, will be resized at use.
     */
}


void ExtendedCCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);      // ue is a POINTER to a UnitigExtension
        ue->setID(i);
        ++i;
    }
    id_init_status = true;
}


void ExtendedCCDBG::print_ids(){
    if (is_id_init() == true){
        std::cout << "[DEBUG] ";
        for (auto &unitig : *this){
            DataAccessor<UnitigExtension>* da = unitig.getData();
            UnitigExtension* ue = da->getData(unitig);
            unsigned id = ue->getID();
            std::cout << id << ", ";
        }
        std::cout << std::endl;
    }
    else
        cerr << "[WARNING] Unitig IDs were not printed because they are not initialized." << endl;
}


/*!
 * \fn      size_t ExtendedCCDBG::count_connected_components()
 * \brief   This function computes the amount of distict connected components in the graph. The connected_components()
 *          function need to be ran successfully to call count_connected_components(). This function is mainly for
 *          debug and test purposes.
 * \return  number of distinct connected components
 */
size_t ExtendedCCDBG::count_connected_components(){
    std::unordered_set<unsigned> unique_set;
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        unique_set.insert(seqan::findSet(UF, ue->getID()));
    }
    return unique_set.size();
}


/*!
 * \fn      bool ExtendedCCDBG::connected_components(const CDBG_Build_opt &graph_options)
 * \brief   This function computes the connected components for at the current state of the graph.
 * \ref     seqan/include/seqan/misc/union_find.h
 * \return  true if successful
 */
bool ExtendedCCDBG::connected_components(const CCDBG_Build_opt &graph_options){
    // initiate UF structure
    if (graph_options.verbose) std::cout << "[VERBOSE] Initiating UNION-FIND" << std::endl;
    resize(UF, (*this).size()+1);   // however, UF needs to be +1 bigger than the graph size. Probably due to start index 1 of IDs
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
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        size_t unitig_id = ue->getID();

        for (auto &it_pre : unitig.getPredecessors()){
            DataAccessor<UnitigExtension>* da = it_pre.getData();
            UnitigExtension* ue = da->getData(it_pre);
            unsigned pre_id = ue->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, pre_id));
        }

        for (auto &it_suc : unitig.getSuccessors()){
            DataAccessor<UnitigExtension>* da = it_suc.getData();
            UnitigExtension* ue = da->getData(it_suc);
            unsigned suc_id = ue->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, suc_id));
        }
    }

    return true;
}


/*!
 * \fn          float ExtendedCCDBG::entropy(const std::string &sequence)
 * \brief       This function computes an entropy for a given string that can be used to filter/mark low complexity
 *              sequences. If all dimers are equaly distributed the entropy is high ("highly chaotic system"), if
 *              all dimers follow a certain pattern the entropy is low ("highly ordered system"). We'd probably like
 *              to mark low entropy unitigs since they have a chance to disrupt/branch the de Bruijn Graph.
 * \remark      Function taken from PopIns.
 * \return      The entropy [0,1] of all binucleotides.
 */
float ExtendedCCDBG::entropy(const std::string &sequence){
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


void ExtendedCCDBG::DFS_cleaner(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        ue->set_undiscovered();
    }
}


void ExtendedCCDBG::DFS_cleaner_seen_only(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        // reset only internal nodes
        if (ue->is_seen())
            ue->set_undiscovered();
    }
}


/*!
 * \fn      inline uint8_t ExtendedCCDBG::whereToGo(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const
 * \brief   This function tests the predecessors P of a unitig u. If the searched unitig src is in P, then
 *          whereToGo() returns GO_FORWARD, denoting the traversal has to continue in the successors of u. If src is
 *          not in P, then whereToGo() returns GO_BACKWARD, denoting the traversal has to continue in the predecessors
 *          of u.
 * \details NOTE: This function is a lowest-level indication where to go, independent of the unitig's orientation. If
 *          the orientation of a unitig is rev-comp, then the result might be GO_BACKWARD while we still consider it
 *          a forward motion with respect to the traversal.
 * \return  DIRECTION
 */
inline uint8_t ExtendedCCDBG::whereToGo(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    uint8_t ret = GO_BACKWARD;
    for (auto &predecessor : um.getPredecessors())
         if (predecessor == src)    /* NOTE: If traversal is within a small-loop, it will always evaluate to true. Therefore loop cases have to be catched during traversal.*/
             ret = GO_FORWARD;
    return ret;
}


/*!
 * \fn      bool ExtendedCCDBG::DFS_Init(const UnitigColorMap< UnitigExtension >& ucm, const bool verbose)
 * \brief   This function initiates the recursion of the directed DFS.
 * \return  bool; 0 for success
 */
bool ExtendedCCDBG::DFS_Init(const UnitigColorMap< UnitigExtension >& ucm, const bool verbose){

    bool ret = 0;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();
    ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

    if (!bw_neighbors.hasPredecessors() && !fw_neighbors.hasSuccessors()){      /* handle singletons */
        // TODO: output sequence of unitig
        ue->set_visited();
    }
    else if(bw_neighbors.hasPredecessors() && fw_neighbors.hasSuccessors()){    /* handle internal nodes */
        // NOTE: just do nothing in this case (internal node)
    }
    else{
        const uint8_t direction = !bw_neighbors.hasPredecessors() ? GO_FORWARD : GO_BACKWARD;

        ue->set_visited();
        if (verbose) cout << "I am setting " << ue->getID() << " to visited." << endl;

        if (direction==GO_BACKWARD){

            for (auto &predecessor : bw_neighbors){

                DataAccessor<UnitigExtension>* da_pre = predecessor.getData();
                UnitigExtension* ue_pre = da_pre->getData(predecessor);

                if (ue_pre->is_undiscovered()){

                    if (verbose) cout << "I am at " << ue->getID() << " and will go backward to " << ue_pre->getID() << endl;

                    //(true == endsHaveSameColors(ucm, predecessor)) ? cout << "SAME colors" << endl : cout << " NOT same colors" << endl;    // TODO delete debug print

                    ret |= DFS_Visit(predecessor, ucm, ucm, verbose);

                    if (verbose) cout << "I jumped back to ID " << ue->getID() << endl;
                }
            }
        }
        else if (direction==GO_FORWARD){

            for (auto &successor : fw_neighbors){

                DataAccessor<UnitigExtension>* da_suc = successor.getData();
                UnitigExtension* ue_suc = da_suc->getData(successor);

                if (ue_suc->is_undiscovered()){

                    if (verbose) cout << "I am at " << ue->getID() << " and will go forward to " << ue_suc->getID() << endl;

                    //(true == endsHaveSameColors(ucm, successor)) ? cout << "SAME colors" << endl : cout << " NOT same colors" << endl;  // TODO delete debug print

                    ret |= DFS_Visit(successor, ucm, ucm, verbose);

                    if (verbose) cout << "I jumped back to ID " << ue->getID() << endl;
                }
            }
        }
        else{
            cerr << "ERROR: " << (direction) << " is not a valid direction" << endl;
            return 1;
        }

        if (verbose) cout << "I am done with " << ue->getID() << endl;
    }

    return ret;
}


/*!
 * \fn      bool ExtendedCCDBG::DFS_Visit(const UnitigColorMap< UnitigExtension >& ucm, const UnitigColorMap< UnitigExtension >& src, const UnitigColorMap< UnitigExtension >& anchor, const bool verbose)
 * \brief   This function executes the recursion of the directed DFS.
 * \return  bool; 0 for success
 */
bool ExtendedCCDBG::DFS_Visit(const UnitigColorMap<UnitigExtension>& ucm,
                              const UnitigColorMap<UnitigExtension>& src,
                              const UnitigColorMap<UnitigExtension>& anchor,
                              const bool verbose){

    bool ret = 0;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    if (whereToGo(ucm, src)==GO_BACKWARD){

        BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();

        // if sink node
        if (!bw_neighbors.hasPredecessors() && !ue->is_visited()){      // TODO: 2nd condition necessary here?
            if (verbose) cout << "I see " << ue->getID() << " has no predecessors and is not visited." << endl;

            DataAccessor<UnitigExtension>* da_anchor = anchor.getData();
            UnitigExtension* ue_anchor = da_anchor->getData(anchor);

            if (verbose) cout << "I will trigger traceback from " << ue->getID() << " to " << ue_anchor->getID() << endl;
            // TODO set sink as visited
            // TODO trigger traceback and return to recursion level above
            return ret;
        }

        ue->set_seen();
        if (verbose) cout << "I am setting " << ue->getID() << " to seen." << endl;

        for (auto &predecessor : bw_neighbors){

            DataAccessor<UnitigExtension>* da_pre = predecessor.getData();
            UnitigExtension* ue_pre = da_pre->getData(predecessor);

            if (ue_pre->is_undiscovered()){

                if (verbose) cout << "I am at " << ue->getID() << " and will go backward to " << ue_pre->getID() << endl;

                ret = DFS_Visit(predecessor, ucm, anchor, verbose);

                if (verbose) cout << "I jumped back to ID " << ue->getID() << endl;
            }
        }
    }

    else{

        ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

        // if sink node
        if (!fw_neighbors.hasSuccessors() && !ue->is_visited()){    // TODO: 2nd condition necessary here?
            if (verbose) cout << "I see " << ue->getID() << " has no successors and is not visited." << endl;

            DataAccessor<UnitigExtension>* da_anchor = anchor.getData();
            UnitigExtension* ue_anchor = da_anchor->getData(anchor);

            if (verbose) cout << "I will trigger traceback from " << ue->getID() << " to " << ue_anchor->getID() << endl;
            // TODO set sink as visited
            // TODO trigger traceback and return to recursion level above
            return ret;
        }

        ue->set_seen();
        if (verbose) cout << "I am setting " << ue->getID() << " to seen." << endl;

        for (auto &successor : fw_neighbors){

            DataAccessor<UnitigExtension>* da_suc = successor.getData();
            UnitigExtension* ue_suc = da_suc->getData(successor);

            if (ue_suc->is_undiscovered()){

                if (verbose) cout << "I am at " << ue->getID() << " and will go forward to " << ue_suc->getID() << endl;

                ret = DFS_Visit(successor, ucm, anchor, verbose);

                if (verbose) cout << "I jumped back to ID " << ue->getID() << endl;
            }
        }
    }

    return ret;
}


bool ExtendedCCDBG::endsHaveSameColors(const UnitigColorMap<UnitigExtension> &observed, const UnitigColorMap<UnitigExtension> &neighbor) const{

    // reverse the return value of whereToGo since this function tells you the opposite direction of where "neighbor" was found as neighbor of observed
    uint8_t direction = (whereToGo(observed, neighbor) == GO_FORWARD) ? GO_BACKWARD : GO_FORWARD;

    size_t nb_colors = this->getNbColors();

    bool same = true;

    if (direction == GO_FORWARD){

        // get kmers
        const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(observed.size - static_cast<size_t>(this->getK()));
        const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(0);

        // get compressed UnitigColor objects
        const UnitigColors* uc_observed_tail = mapping_observed.getData()->getUnitigColors(mapping_observed);
        const UnitigColors* uc_neighbor_head = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);

        for (size_t color_id = 0; color_id != nb_colors; ++color_id){
            same = uc_observed_tail->contains(mapping_observed, color_id) == uc_neighbor_head->contains(mapping_neighbor, color_id);
            if (!same) break;
        }
    }

    else{

        // get kmers
        const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(0);
        const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(neighbor.size - static_cast<size_t>(this->getK()));

        // get compressed UnitigColor objects
        const UnitigColors* uc_observed_head = mapping_observed.getData()->getUnitigColors(mapping_observed);
        const UnitigColors* uc_neighbor_tail = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);

        for (size_t color_id = 0; color_id != nb_colors; ++color_id){
            same = uc_observed_head->contains(mapping_observed, color_id) == uc_neighbor_tail->contains(mapping_neighbor, color_id);
            if (!same) break;
        }
    }

    return same;
}























