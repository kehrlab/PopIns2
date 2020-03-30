#include "ColoredDeBruijnGraph.h"



ExtendedCCDBG::ExtendedCCDBG(int kmer_length, int minimizer_length) :
    ColoredCDBG<UnitigExtension> (kmer_length, minimizer_length),
    id_init_status(false),
    entropy_init_status(false) {
}


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


uint8_t ExtendedCCDBG::traverse(){

    if (!is_id_init()){     // sanity check
        cerr << "[ExtendedCCDBG::traverse] Traversal didn't start because unitig IDs were not initialized." << endl;
        return 0;
    }

    for (auto &ucm : *this){

        if (!is_startnode(ucm)) continue;                       // NOTE: improvement: receive traversal direction from is_startnode() s.t. traverse() doesn't have to find out again

        std::cout << get_unitig_id(ucm) << ": I am a startnode." << std::endl; // DEBUG

        // I think there is no need to color startnodes, since by definition they cannot be part of a cycle.

        if (ucm.getPredecessors().hasPredecessors()){

            std::cout << get_unitig_id(ucm) << ": I jump to my predecessors." << std::endl; // DEBUG

            Traceback tb(VISIT_PREDECESSOR, this->getK());

            if(!DFS(ucm, VISIT_PREDECESSOR, tb)){

                tb.add(ucm.mappedSequenceToString(), true);     // add startnode to final contig

                std::cout << get_unitig_id(ucm) << ": Added start kmer to TB." << std::endl; // DEBUG

            }

            tb.print();
        }
        else{ //ucm.getSuccessors().hasSuccessors()

            std::cout << get_unitig_id(ucm) << ": I jump to my successors." << std::endl; // DEBUG

            Traceback tb(VISIT_SUCCESSOR, this->getK());

            if(!DFS(ucm, VISIT_SUCCESSOR, tb)){

                tb.add(ucm.mappedSequenceToString(), true);     // add startnode to final contig

                std::cout << get_unitig_id(ucm) << ": Added start kmer to TB." << std::endl; // DEBUG

            }

            tb.print();

        }

        reset_dfs_states();

        /*DEBUG*/ std::cout << "" << std::endl;

    }   // end for all unitigs

    return 1;
}


uint8_t ExtendedCCDBG::DFS(const UnitigColorMap<UnitigExtension> &ucm, const direction_t direction, Traceback &tb){

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

                print_unitig_id(ucm);
                print_unitig_seq(ucm);

                std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                tb.add(ucm.mappedSequenceToString());       // add unitig to final contig

                std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                return 0;
            }

            ordered_multimap ranked_predecessors;
            rank_neighbors(ranked_predecessors, ucm, predecessors, VISIT_PREDECESSOR);

            for (auto it = ranked_predecessors.cbegin(); it != ranked_predecessors.cend(); ++it){           // check all ranked neighbors starting from the best color fit

            // traverse further
                for (auto &pre : predecessors){

                    DataAccessor<UnitigExtension>* pre_da = pre.getData();
                    UnitigExtension* pre_data = pre_da->getData(pre);
                    unsigned pre_id = pre_data->getID();

                    unsigned ranked_pre_id = it->second;

                    if(pre_id==ranked_pre_id){

                        if(!DFS(pre, VISIT_PREDECESSOR, tb)){       // if deeper recursion level retuns 0, then stop further traversal here

                            print_unitig_id(ucm);
                            print_unitig_seq(ucm);

                            std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                            tb.add(ucm.mappedSequenceToString());       // add unitig to final contig

                            std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                            return 0;
                        }
                        break;                                  // once the n-th best neighbor was found, (pre_id==ranked_pre_id) won't ever evaluate to true again
                    }
                }
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

                print_unitig_id(ucm);
                print_unitig_seq(ucm);

                std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                tb.add(ucm.mappedSequenceToString());       // add unitig to final contig

                std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                return 0;
            }

            ordered_multimap ranked_successors;
            rank_neighbors(ranked_successors, ucm, successors, VISIT_SUCCESSOR);

            for (auto it = ranked_successors.cbegin(); it != ranked_successors.cend(); ++it){           // check all ranked neighbors starting from the best color fit

            // traverse further
                for (auto &suc : successors){

                    DataAccessor<UnitigExtension>* suc_da = suc.getData();
                    UnitigExtension* suc_data = suc_da->getData(suc);
                    unsigned suc_id = suc_data->getID();

                    unsigned ranked_suc_id = it->second;

                    if(suc_id==ranked_suc_id){

                        if(!DFS(suc, VISIT_SUCCESSOR, tb)){         // if deeper recursion level retuns 0, then stop further traversal here

                            print_unitig_id(ucm);
                            print_unitig_seq(ucm);

                            std::cout << get_unitig_id(ucm) << ": I will jump back." << std::endl; // DEBUG

                            tb.add(ucm.mappedSequenceToString());       // add unitig to final contig

                            std::cout << get_unitig_id(ucm) << ": Added sequence to TB." << std::endl; // DEBUG

                            return 0;
                        }
                        break;                                  // once the n-th best neighbor was found, (suc_id==ranked_suc_id) won't ever evaluate to true again
                    }
                }
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
inline void ExtendedCCDBG::rank_neighbors(ordered_multimap &omm, const UnitigColorMap<UnitigExtension> &ucm, const TNeighbors &neighbors, const direction_t direction) const{

    if (direction==VISIT_PREDECESSOR){

        for (auto &pre : neighbors){

            float overlap = get_neighbor_overlap(ucm, pre);                 // Mind the orientation of the parameter!

            DataAccessor<UnitigExtension>* da = pre.getData();              // NOTE: can I avoid recreating these pointer over and over again?
            UnitigExtension* data = da->getData(pre);

            omm.insert(std::pair<float, unsigned>(overlap, data->getID()));
        }
    }
    else{   // direction==VISIT_SUCCESSOR

        for (auto &suc : neighbors){

            float overlap = get_neighbor_overlap(suc, ucm);                 // Mind the orientation of the parameter!

            DataAccessor<UnitigExtension>* da = suc.getData();              // NOTE: can I avoid recreating these pointer over and over again?
            UnitigExtension* data = da->getData(suc);

            omm.insert(std::pair<float, unsigned>(overlap, data->getID()));
        }
    }
}


inline float ExtendedCCDBG::get_neighbor_overlap(const UnitigColorMap<UnitigExtension> &ucm_to_get_head_from, const UnitigColorMap<UnitigExtension> &ucm_to_get_tail_from) const{
        size_t len = ucm_to_get_tail_from.len;   // I assume this gets me the past-last-kmer index
        //std::cout << len << std::endl;

        const UnitigColorMap<UnitigExtension> k_first = ucm_to_get_head_from.getKmerMapping(0);
        const UnitigColorMap<UnitigExtension> k_last  = ucm_to_get_tail_from.getKmerMapping(len-1);

        const UnitigColors* k_first_colors = k_first.getData()->getUnitigColors(k_first);
        const UnitigColors* k_last_colors  =  k_last.getData()->getUnitigColors(k_last);

        std::vector<bool> k_first_color_bits(this->getNbColors(), false);
        std::vector<bool> k_last_color_bits(this->getNbColors(), false);

        // get color IDs of unitig's first kmer
        UnitigColors::const_iterator cit = k_first_colors->begin(k_first);
        for (; cit != k_first_colors->end(); ++cit)
            k_first_color_bits[cit.getColorID()] = true;

        // get color IDs of unitig's last kmer
        cit = k_last_colors->begin(k_last);
        for (; cit != k_last_colors->end(); ++cit)
            k_last_color_bits[cit.getColorID()] = true;

        // PRINT
        //std::string first_seq = k_first.mappedSequenceToString();
        //std::string last_seq  =  k_last.mappedSequenceToString();
        //std::cout << "START COLORS (" << first_seq << "): " ; print(k_first_color_ids);
        //std::cout << "END COLORS ("   << last_seq  << "): " ; print(k_last_color_ids);

        // calculate Jaccard index
        unsigned numerator = 0;
        unsigned denominator = 0;

        // NOTE: I tried to keep this compliant with the SIMD auto-vectorization of gcc -O3
        // features: support for if-conversion, support for summation reduction
        // unclear: not sure if boolean operations (&&/||) are supported
        size_t i_end = this->getNbColors();
        for (size_t i = 0; i < i_end; ++i){
            numerator   += (k_first_color_bits[i] && k_last_color_bits[i] ? 1 : 0);
            denominator += (k_first_color_bits[i] || k_last_color_bits[i] ? 1 : 0);
        }

        if(!denominator)
            cerr << "[popins2 merge] WARNING: Denominator should never be zero. There has to be at least one color in the graph. Something went wrong!" << endl;

        return (float)numerator / (float)denominator;

        // calculate intersection
        /*
        unsigned count = 0;
        for (unsigned i=0; i < k_first_color_bits.size(); ++i)                  // NOTE: IMPROVEMENT: comparisons can be reduced by using UnitigColors::colorMax(ucm)
            if (k_first_color_bits[i] && k_last_color_bits[i])
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


inline uint8_t ExtendedCCDBG::post_jump_continue_direction(const Kmer &partner, const unsigned lecc_id) const{

    bool go_bw = true;
    bool go_fw = true;

    UnitigMap<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, true> ucm = this->find(partner, true);   // const ucm

    for (auto &pre : ucm.getPredecessors()){

        const DataAccessor<UnitigExtension>* da_pre = pre.getData();
        const UnitigExtension* ue_pre = da_pre->getData(pre);

        if(ue_pre->getLECC() == lecc_id){
            go_bw = false;
            break;
        }
    }

    if (!go_bw)                     // any of the predecessors was associated with the LECC, go on with successors
        return VISIT_SUCCESSOR;

    for (auto &suc : ucm.getSuccessors()){

        const DataAccessor<UnitigExtension>* da_suc = suc.getData();
        const UnitigExtension* ue_suc = da_suc->getData(suc);

        if(ue_suc->getLECC() == lecc_id){
            go_fw = false;
            break;
        }
    }

    // if ( go_bw &&  go_fw) return 0x2 as debug; it means the Kmer has LECC predecessor(s) and successor(s)
    // if (!go_bw && !go_fw) should never happen; Kmer "partner" is always on a LECC border
    return (go_bw && !go_fw) ? VISIT_PREDECESSOR : 0x2;
}






// EOF
