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

        // skip if already LECC-annotated (0 is default LECC identifier)
        if(ue->getLECC() != 0)
            continue;

        // skip high entropy unitigs
        if(ue->getEntropy() >= this->threshold_)
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

    this->lecc_init_status_ = true;

    return LECC_;
}


void LECC_Finder::annotate_recursion(UnitigColorMap<UnitigExtension> &ucm,
                                     const unsigned LECC__){

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    // skip if already LECC-annotated (0 is default LECC identifier), happens in LECC loops
    if(ue->getLECC() != 0)
        return;

    // traversal jumped out of the LECC
    if(ue->getEntropy() >= this->threshold_)
        return;

    // still inside newly discovered LECC
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


inline unsigned LECC_Finder::create_random_color(){
    // RGB color encoding needs a six digit hex number
    // the max six digit hex 0xFFFFFF == 16,777,215 == (2^24)-1 and therefore fits in an unsigned int
    return rand()%16777216;
}


char const * const LECC_Finder::hex_characters = "0123456789ABCDEF";


inline void LECC_Finder::u2hex(std::string &hex, unsigned dec, const unsigned length){
    unsigned pos = length-1;
    while(true){
        unsigned i = dec / 16;
        unsigned r = dec % 16;
        hex[pos] = this->hex_characters[r];
        if(i==0)
            break;
        dec = i;
        --pos;
    }
}


bool LECC_Finder::write(const std::string ofname){
    color_map CM;

    try{
        // FILE
        ofstream ofile(ofname);
        ofile << "ID,Colour" << endl;

        unsigned col;
        std::string col_hex(6, '0');

        srand(time(NULL));
        for (auto &ucm : *g_){
            DataAccessor<UnitigExtension>* da = ucm.getData();
            UnitigExtension* ue = da->getData(ucm);
            unsigned lecc = ue->getLECC();
            unsigned id   = ue->getID();

            // skip if unitig is not in a lecc
            if(lecc==0)
                continue;

            // if lecc has no color assigned yet
            color_map::const_iterator got = CM.find(lecc);
            if(got == CM.end()){
                col = this->create_random_color();
                this->u2hex(col_hex, col);
                ofile << id << ",#" << col_hex << endl;
                // remember color
                CM.insert(color_map_element(lecc,col));
            }
            // if lecc was seen before and therefore already has a color
            else{
                col = got->second;
                this->u2hex(col_hex, col);
                ofile << id << ",#" << col_hex << endl;
            }
        }

        ofile.close();
    }
    catch(std::ofstream::failure &writeErr){
        std::cerr << "\n[ERROR][LECC_Finder::write]: " << writeErr.what() << std::endl;
        return 0;
    }

    return 1;
}


inline bool LECC_Finder::get_borders(border_map_t &border_kmers, const unsigned lecc_id) const{

    for (auto &ucm : *g_){

        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        // only look for unitigs with given LECC ID
        if (ue->getLECC() != lecc_id)
            continue;

        // check if predecessors are borders
        for (auto &pre : ucm.getPredecessors()){

            DataAccessor<UnitigExtension>* da_pre = pre.getData();
            UnitigExtension* ue_pre = da_pre->getData(pre);

            if (ue_pre->getLECC() == 0)
                border_kmers.insert(std::make_pair<Kmer, bool>(pre.getMappedTail(), false));
        }

        // check if successors are borders
        for (auto &suc : ucm.getSuccessors()){

            DataAccessor<UnitigExtension>* da_suc = suc.getData();
            UnitigExtension* ue_suc = da_suc->getData(suc);

            if (ue_suc->getLECC() == 0)
                border_kmers.insert(std::make_pair<Kmer, bool>(suc.getMappedHead(), false));
        }
    }

    if (border_kmers.empty())
        return false;

    return true;
}


inline float LECC_Finder::color_overlap(const Kmer &kmer1, const Kmer &kmer2) const{

    UnitigColorMap<UnitigExtension> ucm_1 = g_->find(kmer1, true);
    UnitigColorMap<UnitigExtension> ucm_2 = g_->find(kmer2, true);

    const UnitigColors* colors_1 = ucm_1.getData()->getUnitigColors(ucm_1);
    const UnitigColors* colors_2 = ucm_2.getData()->getUnitigColors(ucm_2);

    size_t nbColors = g_->getNbColors();

    std::vector<bool> color_bits_1(nbColors, false);
    std::vector<bool> color_bits_2(nbColors, false);

    // get colors IDs and flip the corresponding bits in the color-Bitvector; color IDs are within [0, #color)
    // TODO: this can be done without bitvectors if I compare and update iterators

    UnitigColors::const_iterator cit = colors_1->begin(ucm_1);
    for (; cit != colors_1->end(); ++cit)
        color_bits_1[cit.getColorID()] = true;

    cit = colors_2->begin(ucm_2);   // recycle iterator
    for (; cit != colors_2->end(); ++cit)
        color_bits_2[cit.getColorID()] = true;

    // calculate Jaccard index
    unsigned numerator = 0;
    unsigned denominator = 0;

    // NOTE: I tried to keep this compliant with the SIMD auto-vectorization of gcc -O3
    // features: support for if-conversion, support for summation reduction
    // unclear: not sure if boolean operations (&&/||) are supported
    for (size_t i = 0; i < nbColors; ++i){
        numerator   += (color_bits_1[i] && color_bits_2[i] ? 1 : 0);
        denominator += (color_bits_1[i] || color_bits_2[i] ? 1 : 0);
    }

    if(!denominator)
        cerr << "[popins2 merge] WARNING: Color overlap should never be zero. Something went wrong!" << endl;

    return (float)numerator / (float)denominator;
}


inline void LECC_Finder::reset_dfs_states() const{
    for (auto &ucm : *g_){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        ue->set_undiscovered_fw();
        ue->set_undiscovered_bw();
    }
}


inline void LECC_Finder::check_accessibility(border_map_t &border_kmers, const Kmer &kmer, const unsigned lecc_id) const{

    UnitigColorMap<UnitigExtension> ucm = g_->find(kmer, true);       // unitig of the border Kmer

    for (auto &pre : ucm.getPredecessors()){

        DataAccessor<UnitigExtension>* da_pre = pre.getData();
        UnitigExtension* ue_pre = da_pre->getData(pre);

        if(ue_pre->getLECC() == lecc_id)
            DFS(border_kmers, pre, VISIT_PREDECESSOR);

        reset_dfs_states();
    }

    for (auto &suc : ucm.getSuccessors()){

        DataAccessor<UnitigExtension>* da_suc = suc.getData();
        UnitigExtension* ue_suc = da_suc->getData(suc);

        if(ue_suc->getLECC() == lecc_id)
            DFS(border_kmers, suc, VISIT_SUCCESSOR);

        reset_dfs_states();
    }
}


bool LECC_Finder::DFS(border_map_t &border_kmers, const UnitigColorMap<UnitigExtension> &ucm, const direction_t d) const{

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* data = da->getData(ucm);

    // check if ucm was visited before
    bool is_undiscovered = (d == VISIT_PREDECESSOR ? data->is_undiscovered_bw() : data->is_undiscovered_fw());
    if (!is_undiscovered) return 1;

    // mark current unititg as seen
    d == VISIT_PREDECESSOR ? data->set_seen_bw() : data->set_seen_fw();

    // sink
    if (!data->getLECC()){      // 0 means not in a LECC

        // get border kmer, depending on traversal direction
        border_map_t::iterator got = (d == VISIT_PREDECESSOR ? border_kmers.find(ucm.getMappedTail()) : border_kmers.find(ucm.getMappedHead()));

        // sanity check
        if (got == border_kmers.end()){
            cerr << "[popins2 merge] ERROR: Unitigs immediately outside a LECC should all be member of border_kmers!" << endl;
            return 0;
        }

        // flip accessibility bit
        got->second = true;

        return 1;
    }

    // continue DFS
    if (d == VISIT_PREDECESSOR){
        for (auto &pre : ucm.getPredecessors())
            if(!DFS(border_kmers, pre, VISIT_PREDECESSOR))
                return 0;
    }
    else{   // d == VISIT_SUCCESSOR
        for (auto &suc : ucm.getSuccessors())
            if(!DFS(border_kmers, suc, VISIT_SUCCESSOR))
                return 0;
    }

    return 1;
}


bool LECC_Finder::find_jumps(jump_map_t &jump_map, const unsigned nb_leccs){

    if (!lecc_init_status_){                    // sanity check
        cerr << "[popins2 merge] ERROR: LECC IDs need to be initialized before find_jumps()!" << endl;
        return 0;
    }

    // for every LECC
    for (unsigned i=1; i <= nb_leccs; ++i){

        border_map_t border_kmers;                 // storage for the borders per LECC

        get_borders(border_kmers, i);

        // all every border kmer
        for (auto &border : border_kmers){

            check_accessibility(border_kmers, border.first, i);

            float highscore = 0.0f;
            Kmer jump_target;

            for (auto &border2 : border_kmers){

                if (border2.second == false) continue;   // skip the not accessible partner

                if (border2 == border) continue;         // NOTE this condition has to be evaluated on real data

                float jaccard_index = color_overlap(border.first, border2.first);

                if (jaccard_index > highscore){

                    highscore = jaccard_index;

                    jump_target = border2.first;
                }
            }

            //std::cout << "highscore: " << highscore << '\n';

            if (highscore == 0.0f){
                cerr << "[popins2 merge] WARNING: There was no best jump found for (" << border.first.toString() << ")" << endl;
                continue;
            }

            std::pair<uint64_t,Kmer> jump_pair(border.first.hash(),jump_target);
            jump_map.insert(jump_pair);

            // reset accessibility bits
            for (auto &border : border_kmers)
                border.second = false;

        }   // end kmer border for
    }   // end LECC for

    return 1;
}









// EOF
