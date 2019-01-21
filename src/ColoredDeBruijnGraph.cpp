#include "ColoredDeBruijnGraph.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>


#ifdef DEBUG

void Traceback::printIds() const{
    std::cout << "[";
    for (auto it : this->ids){    
        std::cout << "[";
        for (auto itt : it){
            // output format
            std::cout << itt << ", ";
        }
        std::cout << "]";
    }
    std::cout << "]" << std::endl;
}


void Traceback::printOris() const{
    std::cout << "[";
    for (auto it : this->oris){    
        std::cout << "[";
        for (auto itt : it){
            // output format
            std::cout << (itt ? "+" : "-") << ", ";
        }
        std::cout << "]";
    }
    std::cout << "]" << std::endl;
}


void Traceback::printSeqs() const{
    std::cout << "[";
    for (auto it : this->seqs){    
        std::cout << "[";
        for (auto itt : it){
            // output format
            std::cout << itt << ", ";
        }
        std::cout << "]";
    }
    std::cout << "]" << std::endl;
}


void Traceback::printPathSeqs() const{
    std::cout << "[";
    for (auto it = cbegin(); it != cend(); ++it){
        std::cout << "[";
        for (auto itt : *it){
            // output format
            std::cout << itt << ", ";
        }
        std::cout << "]";
    }
    std::cout << "]" << std::endl;
}

#endif // DEBUG


inline void Traceback::join(const Traceback &t){
    assert(t.seqs.size()==t.ids.size());
#ifdef DEBUG
    for (unsigned i=0; i < t.seqs.size(); ++i){
        this->ids.push_back(t.ids.at(i));
        this->oris.push_back(t.oris.at(i));
        this->seqs.push_back(t.seqs.at(i));
    }
#endif // DEBUG

    for (auto it = t.cbegin(); it != t.cend(); ++it)
        this->push_back(*it);
}


inline void Traceback::cutconcat(string &s, const VSequences &path, const size_t k) const{
    // since the traceback stored the sequences from sink to source, the concatination has to run backwards over the vector
    for (auto rit = path.crbegin(); rit != path.crend(); ++rit){
        if (rit == path.crbegin())
            s+=(*rit);
        else
            s+=rit->substr(k-1);
    }
}


/*!
 * \fn      bool Traceback::write(const std::string &filename, ofstream &ofs) const
 * \return  true if write was successful
 */
bool Traceback::write(ofstream &ofs, const size_t k, size_t &counter) const{
    if (ofs.is_open()){
        // loop through all traceback paths found in one node
        for (auto it=cbegin(); it!=cend(); ++it){
            std::string seq;
            cutconcat(seq, *it, k);
            ++counter;

            ofs << ">" << counter << "\n";
            ofs << seq << "\n";
        }
        return 1;
    }
    else
        cout << "Error: Unable to open file.";
    return 0;
}


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
        std::cout << "[PRINT] ";
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

#ifdef DEBUG
    std::cout << "UF size " << seqan::length(UF._values) << std::endl;
#endif // DEBUG

    // run UF joining
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
 * \fn          inline float ExtendedCCDBG::entropy(const std::string &sequence)
 * \brief       This function computes an entropy for a given string that can be used to filter/mark low complexity
 *              sequences. If all dimers are equaly distributed the entropy is high ("highly chaotic system"), if
 *              all dimers follow a certain pattern the entropy is low ("highly ordered system"). We'd probably like
 *              to mark low entropy unitigs since they have a chance to disrupt/branch the de Bruijn Graph.
 * \remark      Function taken from PopIns.
 * \return      The entropy [0,1] of all binucleotides.
 */
inline float ExtendedCCDBG::entropy(const std::string &sequence){
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


inline void ExtendedCCDBG::DFS_cleaner(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        ue->set_undiscovered_fw();
        ue->set_undiscovered_bw();
    }
}


inline void ExtendedCCDBG::DFS_cleaner_seen_only(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        // reset only internal nodes
        if (ue->is_seen_fw())
            ue->set_undiscovered_fw();
        if (ue->is_seen_bw())
            ue->set_undiscovered_bw();
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
 * \fn      inline uint8_t ExtendedCCDBG::whereFrom(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const
 * \brief   Reverses the direction of whereToGo(). E.g. if the answer is "GO_BACKWARD", then 
 *          this function indicates where to go to reach the source (src).
 *          src  ----------
 *          um          --------->
 * \details This function can also be used to detemine a neighbor's (NBR) orientation with respect to the current unitig (CU) if src=CU and um=NBR.
 * \return  DIRECTION
 */
inline uint8_t ExtendedCCDBG::whereFrom(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    return (whereToGo(um, src) == GO_FORWARD) ? GO_BACKWARD : GO_FORWARD;
}


/*!
 * \fn      Traceback ExtendedCCDBG::DFS_Init(const UnitigColorMap< UnitigExtension >& ucm, const bool verbose)
 * \brief   This function initiates the recursion of the directed DFS.
 * \return  Traceback object
 */
Traceback ExtendedCCDBG::DFS_Init(const UnitigColorMap<UnitigExtension> &ucm, const bool verbose){

    Traceback tb;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);
    if (verbose) cout << "I am starting at " << ue->getID() << "." << endl;

    BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();
    ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

    // --------------------
    //  handle singletons 
    // --------------------
    if (!bw_neighbors.hasPredecessors() && !fw_neighbors.hasSuccessors()){
        // no DFS states to set here

        // Traceback:
#ifdef DEBUG
        Path currentpath;
        currentpath.push_back(ue->getID());
        tb.ids.push_back(currentpath); 

        std::vector<bool> l_oris;
        l_oris.push_back(ucm.strand);
        tb.oris.push_back(l_oris);

        std::vector<std::string> l_seqs;
        l_seqs.push_back(ucm.referenceUnitigToString());
        tb.seqs.push_back(l_seqs);
#endif // DEBUG

        VSequences vseqs;
        ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
        tb.push_back(vseqs);

        tb.recursive_return_status = true;
        return tb;
    }

    // -----------------------
    // handle internal nodes
    // -----------------------
    else if(bw_neighbors.hasPredecessors() && fw_neighbors.hasSuccessors()){
        // NOTE: just do nothing in this case (internal node)
        if (verbose) cout << "Returning from internal node." << endl;
        return tb;
    }
    // ------------------
    // Traverse further
    // ------------------
    else{
        const uint8_t direction = !bw_neighbors.hasPredecessors() ? GO_FORWARD : GO_BACKWARD;

        if (direction==GO_BACKWARD){
            ue->set_seen_bw();
            if (verbose) cout << "I am setting " << ue->getID() << " to seen (bw)." << endl;
            for (auto &predecessor : bw_neighbors){

                DFS_case(ucm, predecessor, tb, verbose);

            }
        }

        else{   // direction==GO_FORWARD
            ue->set_seen_fw();
            if (verbose) cout << "I am setting " << ue->getID() << " to seen (fw)." << endl;
            for (auto &successor : fw_neighbors){

                DFS_case(ucm, successor, tb, verbose);

            }
        }

        ue->set_visited_fw();
        ue->set_visited_bw();
        if (verbose) cout << "I am setting " << ue->getID() << " to visited (both)." << endl;
        if (verbose) cout << "I am done with " << ue->getID() << endl;
    }

    return tb;
}


/*!
 * \fn      Traceback ExtendedCCDBG::DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
 *                                             const UnitigColorMap<UnitigExtension> &src,
 *                                             const uint8_t src_direction,
 *                                             const bool verbose)
 * \brief   This function executes the recursion of the directed DFS.
 * \param   ucm is the current node during traversal
 * \param   src_direction is the direction you'd need to go if you wanted to go back from ucm to src. Therefore, for further traversal
 *          you need to follow the opposite direction.
 * \return  Traceback object
 */
Traceback ExtendedCCDBG::DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
                                   const UnitigColorMap<UnitigExtension> &src,
                                   const uint8_t src_direction,
                                   const bool verbose){
    Traceback tb;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    uint8_t traversal_direction = (src_direction==GO_BACKWARD) ? GO_FORWARD : GO_BACKWARD;  // NOTE: I guess this could be avoided if I'd pass traversal direction directly in DFS_visit

    // -------------------------
    // |  if traverse forward  |
    // -------------------------
    if (traversal_direction==GO_BACKWARD){
        if (verbose) cout << "I am setting " << ue->getID() << " to seen (bw)." << endl;
        ue->set_seen_bw();  // NOTE: I keep this marked since we (in it's current version) only report one path per source-sink note

        BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();

        // ------------------
        // |  if sink node  |
        // ------------------
        if ((!bw_neighbors.hasPredecessors() && !ue->is_visited_bw()) || !endsHaveCommonColor(src, ucm, verbose)){   // visited check could be bw/fw; visited check necessary to avoid RC path
            if (verbose) cout << "I see " << ue->getID() << " has no predecessors and is not visited." << endl;
            if (verbose) cout << "I will trigger traceback from " << ue->getID() << endl;

            // Traceback:
#ifdef DEBUG
            Path currentpath;
            currentpath.push_back(ue->getID());
            tb.ids.push_back(currentpath);

            std::vector<bool> l_oris;
            l_oris.push_back(ucm.strand);
            tb.oris.push_back(l_oris);

            std::vector<std::string> l_seqs;
            l_seqs.push_back(ucm.referenceUnitigToString());
            tb.seqs.push_back(l_seqs);
#endif // DEBUG

            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);

            tb.recursive_return_status = true;
            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------
        for (auto &predecessor : bw_neighbors){

            DFS_case(ucm, predecessor, tb, verbose);

        }
    }

    // --------------------------
    // |  if traverse backward  |
    // --------------------------
    else{   // traversal_direction==GO_FORWARD
        if (verbose) cout << "I am setting " << ue->getID() << " to seen (fw)." << endl;
        ue->set_seen_fw();  // NOTE: I keep this marked since we (in it's current version) only report one path per source-sink note

        ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

        // ------------------
        // |  if sink node  |
        // ------------------
        if ((!fw_neighbors.hasSuccessors() && !ue->is_visited_fw()) || !endsHaveCommonColor(src, ucm, verbose)){   // visited check could be bw/fw; visited check necessary to avoid RC path
            if (verbose) cout << "I see " << ue->getID() << " has no successors and is not visited." << endl;
            if (verbose) cout << "I will trigger traceback from " << ue->getID() << endl;

            // Traceback:
#ifdef DEBUG
            Path currentpath;
            currentpath.push_back(ue->getID());
            tb.ids.push_back(currentpath);

            std::vector<bool> l_oris;
            l_oris.push_back(ucm.strand);
            tb.oris.push_back(l_oris);

            std::vector<std::string> l_seqs;
            l_seqs.push_back(ucm.referenceUnitigToString());
            tb.seqs.push_back(l_seqs);
#endif // DEBUG

            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);

            tb.recursive_return_status = true;
            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------
        for (auto &successor : fw_neighbors){

            DFS_case(ucm, successor, tb, verbose);

        }
    }

    return tb;
}


/*!
 * \fn      void ExtendedCCDBG::DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
 *                                       const UnitigColorMap<UnitigExtension> &neighbor,
 *                                       Traceback &tb,
 *                                       const bool verbose)
 * \brief   This function contains the code to determine the neighbor's (pre/suc) orientation, the recursion call of the according case and
 *          the management of a (recursively) returned traceback instance.
 * \return  void
 */
void ExtendedCCDBG::DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
                             const UnitigColorMap<UnitigExtension> &neighbor,
                             Traceback &tb,
                             const bool verbose){

#ifdef DEBUG
    DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
    UnitigExtension* ucm_ue = ucm_da->getData(ucm);
#endif // DEBUG

    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);

    if (whereFrom(neighbor, ucm)==GO_BACKWARD){
        /*  Case 3:                             Case 1:
        *   SRC ------->                OR              -------> SRC
        *           -------> SUC                PRE <--------
        */
        if (neighbor_ue->is_undiscovered_fw()){
            // Recursion
#ifdef DEBUG
            if (verbose) cout << "I am at " << ucm_ue->getID() << " and will go forward to " << neighbor_ue->getID() << endl;
#endif // DEBUG
            Traceback returned_tb = DFS_Visit(neighbor, ucm, GO_BACKWARD, verbose);
#ifdef DEBUG
            if (verbose) cout << "I jumped back to ID " << ucm_ue->getID() << endl;
#endif // DEBUG

            // Traceback
            if (returned_tb.recursive_return_status){
#ifdef DEBUG
                for (unsigned i=0; i < returned_tb.seqs.size(); ++i){
                    returned_tb.seqs.at(i).push_back(ucm.referenceUnitigToString());
                    returned_tb.oris.at(i).push_back(ucm.strand);
                    returned_tb.ids.at(i).push_back(ucm_ue->getID());
                }
                assert(returned_tb.empty() == false)
#endif // DEBUG
                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);
            }
        }
    }
    else{   // whereFrom(neighbor, ucm)==GO_FORWARD
        /*  Case 4:                             Case 2:
        *   SRC ------->                OR              -------> SRC
        *           <------- SUC                PRE ------->
        */
        if (neighbor_ue->is_undiscovered_bw()){
            // Recursion
#ifdef DEBUG
            if (verbose) cout << "I am at " << ucm_ue->getID() << " and will go forward to " << neighbor_ue->getID() << endl;
#endif // DEBUG
            Traceback returned_tb = DFS_Visit(neighbor, ucm, GO_FORWARD, verbose);
#ifdef DEBUG
            if (verbose) cout << "I jumped back to ID " << ucm_ue->getID() << endl;
#endif // DEBUG

            // Traceback recursion
            if (returned_tb.recursive_return_status){
#ifdef DEBUG
                for (unsigned i=0; i < returned_tb.seqs.size(); ++i){
                    returned_tb.seqs.at(i).push_back(ucm.referenceUnitigToString());
                    returned_tb.oris.at(i).push_back(ucm.strand);
                    returned_tb.ids.at(i).push_back(ucm_ue->getID());
                }
                assert(returned_tb.empty() == false)
#endif // DEBUG
                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);
            }
        }
    }
}


/*!
 * \fn      inline bool ExtendedCCDBG::endsHaveCommonColor(const UnitigColorMap<UnitigExtension> &observed,
 *                                                         const UnitigColorMap<UnitigExtension> &neighbor) const
 * \details This function compares the last kmer's color vector of two consecutive unitigs (observed OBS, neighbor NBR), i.e.
 *          OBS -------        -------- NBR
 *                    X               A
 *                    Y  --->         B         if (X==A) || (Y==B) || (Z==C) then return true
 *                    Z               C
 * \return  bool
 */
inline bool ExtendedCCDBG::endsHaveCommonColor(const UnitigColorMap<UnitigExtension> &observed,
                                               const UnitigColorMap<UnitigExtension> &neighbor, 
                                               const bool verbose) const{
    bool rValue = false;
    size_t nb_colors = this->getNbColors();

    // get unitig orientations
    uint8_t observed_to_neighbor_direction = whereFrom(observed, neighbor);
    uint8_t neighbor_to_observed_direction = whereFrom(neighbor, observed);

    if (observed_to_neighbor_direction == GO_FORWARD){
        if(neighbor_to_observed_direction == GO_FORWARD){
            // CASE 1:      OBS -----> <----- NBR
            // Comp. :               X      X
            // get kmers
            const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(observed.size - static_cast<size_t>(this->getK()));
            const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(0);
            // get compressed UnitigColor objects
            const UnitigColors* uc_observed_tail = mapping_observed.getData()->getUnitigColors(mapping_observed);
            const UnitigColors* uc_neighbor_head = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);

            for (size_t color_id = 0; color_id != nb_colors; ++color_id){
                rValue = uc_observed_tail->contains(mapping_observed, color_id) == uc_neighbor_head->contains(mapping_neighbor, color_id);
                if (rValue) break;
            }
        }

        if(neighbor_to_observed_direction == GO_BACKWARD){
            // CASE 2:      OBS -----> -----> NBR
            // Comp. :               X      X
            // get kmers
            const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(observed.size - static_cast<size_t>(this->getK()));
            const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(neighbor.size - static_cast<size_t>(this->getK()));
            // get compressed UnitigColor objects
            const UnitigColors* uc_observed_tail = mapping_observed.getData()->getUnitigColors(mapping_observed);
            const UnitigColors* uc_neighbor_tail = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);
 
            for (size_t color_id = 0; color_id != nb_colors; ++color_id){
                rValue = uc_observed_tail->contains(mapping_observed, color_id) == uc_neighbor_tail->contains(mapping_neighbor, color_id);
                if (rValue) break;
            }
        }
    }

    else{   // observed_to_neighbor_direction == GO_BACKWARD
        if(neighbor_to_observed_direction == GO_FORWARD){
            // CASE 3:      NBR -----> -----> OBS
            // Comp. :          X      X
            // get kmers
            const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(0);
            const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(0);
            // get compressed UnitigColor objects
            const UnitigColors* uc_observed_head = mapping_observed.getData()->getUnitigColors(mapping_observed);
            const UnitigColors* uc_neighbor_head = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);

            for (size_t color_id = 0; color_id != nb_colors; ++color_id){
                rValue = uc_observed_head->contains(mapping_observed, color_id) == uc_neighbor_head->contains(mapping_neighbor, color_id);
                if (rValue) break;
            }
        }

        if(neighbor_to_observed_direction == GO_BACKWARD){
            // CASE 4:      NBR <----- -----> OBS
            // Comp. :          X      X
            // get kmers
            const const_UnitigColorMap<UnitigExtension> mapping_observed = observed.getKmerMapping(0);
            const const_UnitigColorMap<UnitigExtension> mapping_neighbor = neighbor.getKmerMapping(neighbor.size - static_cast<size_t>(this->getK()));
            // get compressed UnitigColor objects
            const UnitigColors* uc_observed_head = mapping_observed.getData()->getUnitigColors(mapping_observed);
            const UnitigColors* uc_neighbor_tail = mapping_neighbor.getData()->getUnitigColors(mapping_neighbor);

            for (size_t color_id = 0; color_id != nb_colors; ++color_id){
                rValue = uc_observed_head->contains(mapping_observed, color_id) == uc_neighbor_tail->contains(mapping_neighbor, color_id);
                if (rValue) break;
            }
        }
    }

    if (verbose)
        if (!rValue) cout << "Color-criterion stopped traversal." << endl;

    return rValue;
}


/*!
 * \fn      bool ExtendedCCDBG::merge(const CCDBG_Build_opt &opt)
 * \return  bool; true if successful
 */
bool ExtendedCCDBG::merge(const CCDBG_Build_opt &opt){
    /* sanity check */
    if (!this->is_id_init())
        return false;
    
    Traceback tb;

    size_t sv_counter = 0;
    std::string sv_filename = "contigs.fa";
    
    ofstream ofs(sv_filename, std::ofstream::out);

    if (!ofs.is_open()){
        cerr << "Error: Couldn't open ofstream for contig file." << endl;
        return false;
    }

    for (auto &unitig : *this){
        tb = DFS_Init(unitig, opt.verbose);
        if (tb.recursive_return_status){
#ifdef DEBUG
            tb.printIds();
            tb.printOris();
            tb.printSeqs();
#endif // DEBUG

            if (!tb.write(ofs, opt.k, sv_counter))
                return false;
        }
        DFS_cleaner_seen_only();
    }

    ofs.close();

    return true;
}











