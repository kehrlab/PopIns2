#include "ColoredDeBruijnGraph.h"



void ExtendedCCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);      // ue is a POINTER to a UnitigExtension
        ue->setID(i);
        ++i;
    }
    this->id_init_status = true;
}


void ExtendedCCDBG::print_ids(){
    if (!is_id_init()){
        cerr << "[WARNING] Unitig IDs were not printed because they are not initialized." << endl;
        return;
    }

    std::cout << "[PRINT] ";
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        unsigned id = ue->getID();
        std::cout << id << ", ";
    }
    std::cout << std::endl;
}


uint8_t ExtendedCCDBG::traverse(const int setcover_threshold){

    // sanity checks
    if (!is_id_init()){
        cerr << "[ExtendedCCDBG::traverse] Traversal didn't start because unitig IDs were not initialized." << endl;
        return 0;
    }

    if (!entropy_init_status){
        cerr << "[ExtendedCCDBG::traverse] Traversal didn't start because unitig IDs were not annotated with entropy." << endl;
        return 0;
    }

    if (_jump_map_ptr==NULL){
        cerr << "[ExtendedCCDBG::traverse] Traversal didn't start because graph got no jump_map_t* assigned." << endl;
        return 0;
    }

    Setcover sc(setcover_threshold);

    // main routine
    for (auto &ucm : *this){

        if (!is_startnode(ucm)) continue;                       // NOTE: improvement: receive traversal direction from is_startnode() s.t. traverse() doesn't have to find out again

        std::cout << get_unitig_id(ucm) << ": I am a startnode." << std::endl; // DEBUG

        // I think there is no need to color startnodes, since by definition they cannot be part of a cycle.

        if (ucm.getPredecessors().hasPredecessors()){

            std::cout << get_unitig_id(ucm) << ": I jump to my predecessors." << std::endl; // DEBUG

            Traceback tb(VISIT_PREDECESSOR, this->getK());

            if(!DFS(ucm, VISIT_PREDECESSOR, tb, sc)){

                tb.add(ucm.referenceUnitigToString(), true);     // add startnode to final contig

                // don't sc.add(ucm) here because the current ucm is added in ExtendedCCDBG::DFS()

                std::cout << get_unitig_id(ucm) << ": Added start kmer to TB." << std::endl; // DEBUG
            }

            if(sc.test())
                tb.print();
        }

        else{ //ucm.getSuccessors().hasSuccessors()

            std::cout << get_unitig_id(ucm) << ": I jump to my successors." << std::endl; // DEBUG

            Traceback tb(VISIT_SUCCESSOR, this->getK());

            if(!DFS(ucm, VISIT_SUCCESSOR, tb, sc)){

                tb.add(ucm.referenceUnitigToString(), true);     // add startnode to final contig

                // don't sc.add(ucm) here because the current ucm is added in ExtendedCCDBG::DFS()

                std::cout << get_unitig_id(ucm) << ": Added start kmer to TB." << std::endl; // DEBUG
            }

            if(sc.test())
                tb.print();
        }

        reset_dfs_states();

        /*DEBUG*/ std::cout << "" << std::endl;

    }   // end for all unitigs

    return 1;
}


uint8_t ExtendedCCDBG::DFS(const UnitigColorMap<UnitigExtension> &ucm, const direction_t direction, Traceback &tb, Setcover &sc, const bool jumped){

    std::cout << get_unitig_id(ucm) << ": I am a neighbor." << std::endl; // DEBUG

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* data = da->getData(ucm);

    if(direction==VISIT_PREDECESSOR){

        if(data->is_undiscovered_bw()){

            data->set_seen_bw();

            std::cout << get_unitig_id(ucm) << ": [PRE] Now I am the current state." << std::endl; // DEBUG

            BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> predecessors = ucm.getPredecessors();

            // sink node
            if (!predecessors.hasPredecessors()){

                std::cout << get_unitig_id(ucm) << ": I am a sink." << std::endl; // DEBUG
                std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                jumped ? tb.addFullSink(ucm.referenceUnitigToString()) : tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                sc.add(ucm);

                std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                return 0;
            }

            ordered_multimap_t ranked_predecessors;
            rank_neighbors(ranked_predecessors, ucm, predecessors, VISIT_PREDECESSOR);

            /*DEBUG
            std::cout << "*** ranked predecessors ***" << std::endl;
            for (auto &t : ranked_predecessors)
                cout << t.first << " => " << t.second << " (float => partner ID)" << endl;
            DEBUG*/

            for (auto it = ranked_predecessors.cbegin(); it != ranked_predecessors.cend(); ++it){           // check all ranked neighbors starting from the best color fit

                bool _was_direct_neighbor = false;

                // traverse predecessors further (direct neighbors)
                for (auto &pre : predecessors){

                    DataAccessor<UnitigExtension>* pre_da = pre.getData();
                    UnitigExtension* pre_data = pre_da->getData(pre);
                    unsigned pre_id = pre_data->getID();

                    unsigned ranked_pre_id = it->second;

                    if(pre_id==ranked_pre_id){

                        _was_direct_neighbor = true;

                        if(!DFS(pre, VISIT_PREDECESSOR, tb, sc)){       // if deeper recursion level retuns 0, then stop further traversal here

                            std::cout << get_unitig_id(ucm) << ": I will step back." << std::endl; // DEBUG

                            tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                            sc.add(ucm);

                            std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                            return 0;
                        }
                        break;                                  // once the n-th best neighbor was found, (pre_id==ranked_pre_id) won't ever evaluate to true again
                    }
                }

                // traverse predecessors further (jump)
                if (!_was_direct_neighbor){

                    uint64_t _border_hash = ucm.getMappedHead().rep().hash();

                    jump_map_t::const_iterator cit = _jump_map_ptr->find(_border_hash);

                    if (cit == _jump_map_ptr->cend()){
                        cerr << "[popins2 merge] WARNING: ExtendedCCDBG::DFS() couldn't find a kmer to jump to." << endl;
                        return 1;       // Look for another way then. TODO: return better error codes in DFS(), 1 is not ideal here.
                    }

                    const Kmer partner_kmer = cit->second;

                    const UnitigColorMap<UnitigExtension> partner_unitig = this->find(partner_kmer, true);

                    std::cout << get_unitig_id(ucm) << ": I will jump over a LECC to " << get_unitig_id(partner_unitig) << std::endl; // DEBUG

                    // TODO: catch error in post_jump_continue_direction() here

                    // TEST: is post_jump_continue_direction() necessary?
                    if(!DFS(partner_unitig, post_jump_continue_direction(partner_unitig), tb, sc, true)){       // if deeper recursion level retuns 0, then stop further traversal here

                        std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                        tb.addN();

                        std::cout << get_unitig_id(ucm) << ": Added Ns to TB." << std::endl; // DEBUG

                        tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                        sc.add(ucm);

                        std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                        return 0;
                    }
                }

                // if control flow reaches this point, for-loop continues with next best ranked predecessor
            }
        }
    }
    else{   // if direction==VISIT_SUCCESSOR

        if(data->is_undiscovered_fw()){

            data->set_seen_fw();

            std::cout << get_unitig_id(ucm) << ": [SUC] Now I am the current state." << std::endl; // DEBUG

            ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> successors = ucm.getSuccessors();

            // sink node
            if (!successors.hasSuccessors()){

                std::cout << get_unitig_id(ucm) << ": I am a sink." << std::endl; // DEBUG
                std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                jumped ? tb.addFullSink(ucm.referenceUnitigToString()) : tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                sc.add(ucm);

                std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                return 0;
            }

            ordered_multimap_t ranked_successors;
            rank_neighbors(ranked_successors, ucm, successors, VISIT_SUCCESSOR);

            /*DEBUG
            std::cout << "*** ranked successors ***" << std::endl;
            for (auto &t : ranked_successors)
                cout << t.first << " => " << t.second << " (float => partner ID)" << endl;
            DEBUG*/

            for (auto it = ranked_successors.cbegin(); it != ranked_successors.cend(); ++it){           // check all ranked neighbors starting from the best color fit

                bool _was_direct_neighbor = false;

                // traverse successors further (direct neighbors)
                for (auto &suc : successors){

                    DataAccessor<UnitigExtension>* suc_da = suc.getData();
                    UnitigExtension* suc_data = suc_da->getData(suc);
                    unsigned suc_id = suc_data->getID();

                    unsigned ranked_suc_id = it->second;

                    if(suc_id==ranked_suc_id){

                        _was_direct_neighbor = true;

                        if(!DFS(suc, VISIT_SUCCESSOR, tb, sc)){         // if deeper recursion level retuns 0, then stop further traversal here

                            std::cout << get_unitig_id(ucm) << ": I will step back." << std::endl; // DEBUG

                            tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                            sc.add(ucm);

                            std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                            return 0;
                        }
                        break;                                  // once the n-th best neighbor was found, (suc_id==ranked_suc_id) won't ever evaluate to true again
                    }
                }

                // traverse successors further (jump)
                if(!_was_direct_neighbor){
                    uint64_t _border_hash = ucm.getMappedTail().rep().hash();

                    jump_map_t::const_iterator cit = _jump_map_ptr->find(_border_hash);

                    if (cit == _jump_map_ptr->cend()){
                        cerr << "[popins2 merge] WARNING: ExtendedCCDBG::DFS() couldn't find a kmer to jump to." << endl;
                        return 1;       // Look for another way then. TODO: return better error codes in DFS(), 1 is not ideal here.
                    }

                    const Kmer partner_kmer = cit->second;

                    const UnitigColorMap<UnitigExtension> partner_unitig = this->find(partner_kmer, true);

                    std::cout << get_unitig_id(ucm) << ": I will jump over a LECC to " << get_unitig_id(partner_unitig) << std::endl; // DEBUG

                    // TODO: catch error in post_jump_continue_direction() here

                    // TEST: is post_jump_continue_direction() necessary?
                    if(!DFS(partner_unitig, post_jump_continue_direction(partner_unitig), tb, sc, true)){       // if deeper recursion level retuns 0, then stop further traversal here

                        std::cout << get_unitig_id(ucm) << ": I will step back." << std::endl; // DEBUG

                        tb.addN();

                        std::cout << get_unitig_id(ucm) << ": Added Ns to TB." << std::endl; // DEBUG

                        tb.add(ucm.referenceUnitigToString());       // add unitig to final contig

                        sc.add(ucm);

                        std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                        return 0;
                    }
                }

                // if control flow reaches this point, for-loop continues with next best ranked successor
            }
        }
    }

    return 1;                                           // let higher recursion level keep on searching (just jump back)
}


inline void ExtendedCCDBG::reset_dfs_states(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        ue->set_undiscovered_fw();
        ue->set_undiscovered_bw();
    }
}


inline bool ExtendedCCDBG::is_startnode(const UnitigColorMap<UnitigExtension> &ucm){
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    return
        ucm.len>2                                                                                       // be longer than 2 kmers and
        &&
        (( ucm.getPredecessors().hasPredecessors() && !ucm.getSuccessors().hasSuccessors() ) ||         // have only predecessors or
         (!ucm.getPredecessors().hasPredecessors() &&  ucm.getSuccessors().hasSuccessors() ))           // have only successors
        &&
        ue->getLECC()==0u;                                                                              // is not part of an LECC
}


template <typename TNeighbors>
inline void ExtendedCCDBG::rank_neighbors(ordered_multimap_t &omm, const UnitigColorMap<UnitigExtension> &ucm, const TNeighbors &neighbors, const direction_t direction){

    DataAccessor<UnitigExtension>* da;
    UnitigExtension* data;

    float overlap;

    if (direction==VISIT_PREDECESSOR){

        for (auto &pre : neighbors){

            da = pre.getData();
            data = da->getData(pre);

            // if neighbor is in LECC consider a jump
            if (data->getLECC()){

                const uint64_t border_hash = ucm.getMappedHead().rep().hash();      // ucm is border here

                /*DEBUG
                cout << "----------------------------"                              << endl;
                cout << "BORDER   ID: " << get_unitig_id(ucm)                       << endl;
                cout << "BORDER KMER: " << ucm.getMappedHead().rep().toString()     << endl;
                cout << "BORDER HASH: " << border_hash                              << endl;
                cout << "----------------------------"                              << endl;
                DEBUG*/

                // assume having a pointer jump_map_ptr to a jump_map_t (std::unordered_map<uint64_t, Kmer>) in this scope
                std::unordered_map<uint64_t, Kmer>::const_iterator cit = _jump_map_ptr->find(border_hash);

                if (cit == _jump_map_ptr->cend()){
                    cerr << "[popins2 merge] WARNING: ExtendedCCDBG::rank_neighbors() couldn't find a kmer to jump to." << endl;
                    continue;           // don't consider the "pre" inside the LECC any further
                }

                const Kmer partner_kmer = cit->second;

                const UnitigColorMap<UnitigExtension> partner_unitig = this->find(partner_kmer, true);

                overlap = get_neighbor_overlap(ucm, partner_unitig);

                da = partner_unitig.getData();
                data = da->getData(partner_unitig);

                omm.insert(std::pair<float, unsigned>(overlap, data->getID()));

                continue;
            }

            overlap = get_neighbor_overlap(ucm, pre);                 // Mind the orientation of the parameter!

            omm.insert(std::pair<float, unsigned>(overlap, data->getID()));
        }
    }
    else{   // direction==VISIT_SUCCESSOR

        for (auto &suc : neighbors){

            da = suc.getData();
            data = da->getData(suc);

            // if neighbor is in LECC consider a jump
            if (data->getLECC()){

                const uint64_t border_hash = ucm.getMappedTail().rep().hash();      // ucm is border here

                /*DEBUG
                cout << "----------------------------"                              << endl;
                cout << "BORDER   ID: " << get_unitig_id(ucm)                       << endl;
                cout << "BORDER KMER: " << ucm.getMappedHead().rep().toString()     << endl;
                cout << "BORDER HASH: " << border_hash                              << endl;
                cout << "----------------------------"                              << endl;
                DEBUG*/

                // assume having a pointer jump_map_ptr to a jump_map_t (std::unordered_map<uint64_t, Kmer>) in this scope
                std::unordered_map<uint64_t, Kmer>::const_iterator cit = _jump_map_ptr->find(border_hash);

                if (cit == _jump_map_ptr->cend()){
                    cerr << "[popins2 merge] WARNING: ExtendedCCDBG::rank_neighbors() couldn't find a kmer to jump to." << endl;
                    continue;           // don't consider the "suc" inside the LECC any further
                }

                const Kmer partner_kmer = cit->second;

                const UnitigColorMap<UnitigExtension> partner_unitig = this->find(partner_kmer, true);

                overlap = get_neighbor_overlap(partner_unitig, ucm);

                da = partner_unitig.getData();
                data = da->getData(partner_unitig);

                omm.insert(std::pair<float, unsigned>(overlap, data->getID()));

                continue;
            }

            overlap = get_neighbor_overlap(suc, ucm);                 // Mind the orientation of the parameter!

            omm.insert(std::pair<float, unsigned>(overlap, data->getID()));
        }
    }
}


inline float ExtendedCCDBG::get_neighbor_overlap(const UnitigColorMap<UnitigExtension> &extract_head, const UnitigColorMap<UnitigExtension> &extract_tail) {

        // --------------------------------------------     // NOTE: this part should be done better

        const Kmer head = extract_head.getMappedHead().rep();
        const Kmer tail = extract_tail.getMappedTail().rep();

        const UnitigColorMap<UnitigExtension> head_ucm = this->find(head, true);
        const UnitigColorMap<UnitigExtension> tail_ucm = this->find(tail, true);

        // --------------------------------------------

        const UnitigColors* head_colors = head_ucm.getData()->getUnitigColors(head_ucm);
        const UnitigColors* tail_colors = tail_ucm.getData()->getUnitigColors(tail_ucm);

        std::vector<bool> head_color_bits(this->getNbColors(), false);
        std::vector<bool> tail_color_bits(this->getNbColors(), false);

        // get color IDs of unitig's first kmer
        UnitigColors::const_iterator cit = head_colors->begin(head_ucm);
        for (; cit != head_colors->end(); ++cit)
            head_color_bits[cit.getColorID()] = true;

        // get color IDs of unitig's last kmer
        cit = tail_colors->begin(tail_ucm);
        for (; cit != tail_colors->end(); ++cit)
            tail_color_bits[cit.getColorID()] = true;

        /*DEBUG
        std::cout << "HEAD COLORS (OF UNITIG ID " << get_unitig_id(head_ucm) << "): " ; prettyprint::print(head_color_bits);
        std::cout << "TAIL COLORS (OF UNITIG ID " << get_unitig_id(tail_ucm) << "): " ; prettyprint::print(tail_color_bits);
        DEBUG*/

        // calculate Jaccard index
        unsigned numerator = 0;
        unsigned denominator = 0;

        // NOTE: I tried to keep this compliant with the SIMD auto-vectorization of gcc -O3
        // features: support for if-conversion, support for summation reduction
        // unclear: not sure if boolean operations (&&/||) are supported
        size_t i_end = this->getNbColors();
        for (size_t i = 0; i < i_end; ++i){
            numerator   += (head_color_bits[i] && tail_color_bits[i] ? 1 : 0);
            denominator += (head_color_bits[i] || tail_color_bits[i] ? 1 : 0);
        }

        if(!denominator)
            cerr << "[popins2 merge] WARNING: Denominator should never be zero. There has to be at least one color in the graph." << endl;

        return (float)numerator / (float)denominator;

        // calculate intersection
        /*
        unsigned count = 0;
        for (unsigned i=0; i < head_color_bits.size(); ++i)                  // NOTE: IMPROVEMENT: comparisons can be reduced by using UnitigColors::colorMax(ucm)
            if (head_color_bits[i] && tail_color_bits[i])
                ++count;

        return count;
        */
}


void ExtendedCCDBG::init_entropy(){

    for (auto &unitig : *this){

        const std::string unitig_s = unitig.referenceUnitigToString();

        const float entropy = this->entropy(unitig_s);

        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        ue->setEntropy(entropy);
    }

    this->entropy_init_status = true;
}


inline float ExtendedCCDBG::entropy(const std::string &sequence){
    // create a dictionary counting the occurrence of all dinucleotides
    std::unordered_map<std::string, unsigned> diCounts(16);

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
    float entropy = 0;

    for(std::unordered_map<std::string,unsigned>::const_iterator it = diCounts.cbegin(); it != diCounts.cend(); ++it){

        if (it->second == 0) continue;

        float p = float(it->second) / counted;

        entropy -= p * seqan::log(p) / seqan::log(2);
    }

    return entropy / 4;
}


inline uint8_t ExtendedCCDBG::post_jump_continue_direction(const UnitigColorMap<UnitigExtension> &ucm) const{

    bool go_bw = true;
    bool go_fw = true;

    for (auto &pre : ucm.getPredecessors()){

        const DataAccessor<UnitigExtension>* da_pre = pre.getData();
        const UnitigExtension* ue_pre = da_pre->getData(pre);

        if(ue_pre->getLECC() != 0){
            go_bw = false;
            break;
        }
    }

    if (!go_bw)                     // any of the predecessors was associated with the LECC, go on with successors
        return VISIT_SUCCESSOR;

    for (auto &suc : ucm.getSuccessors()){

        const DataAccessor<UnitigExtension>* da_suc = suc.getData();
        const UnitigExtension* ue_suc = da_suc->getData(suc);

        if(ue_suc->getLECC() != 0){
            go_fw = false;
            break;
        }
    }

    // if ( go_bw &&  go_fw) return 0x2 as debug; it means the Kmer has LECC predecessor(s) and successor(s)
    // if (!go_bw && !go_fw) should never happen; Kmer "partner" is always on a LECC border
    return (go_bw && !go_fw) ? VISIT_PREDECESSOR : 0x2;
}



// EOF
