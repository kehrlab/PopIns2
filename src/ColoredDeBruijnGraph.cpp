#include "ColoredDeBruijnGraph.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <map>



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


size_t ExtendedCCDBG::count_connected_components(){
    std::unordered_set<unsigned> unique_set;
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        unique_set.insert(seqan::findSet(UF, ue->getID()));
    }
    return unique_set.size();
}


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


inline uint8_t ExtendedCCDBG::whereToGo(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    uint8_t ret = GO_BACKWARD;
    for (auto &predecessor : um.getPredecessors())
         if (predecessor == src)    /* NOTE: If traversal is within a small-loop, it will always evaluate to true. Therefore loop cases have to be catched during traversal.*/
             ret = GO_FORWARD;
    return ret;
}


inline uint8_t ExtendedCCDBG::whereFrom(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    return (whereToGo(um, src) == GO_FORWARD) ? GO_BACKWARD : GO_FORWARD;
}


bool ExtendedCCDBG::merge(const CCDBG_Build_opt &opt){
    // SANITY CHECK(S)
    if (!this->is_id_init())
        return false;

    // I/O
    std::string sv_filename = "contigs.fa";
    ofstream ofs(sv_filename, std::ofstream::out);
    if (!ofs.is_open()){cerr << "Error: Couldn't open ofstream for contig file." << endl; return false;}

    // DFS
    Setcover<> sc;
    size_t sv_counter = 0;
    for (auto &unitig : *this){
            if (opt.verbose) cout << " -------------------------------- " << endl;
            Traceback tb = DFS_Init_bidirectional(unitig, opt.verbose, sc);
            if (tb.recursive_return_status)
                if (!tb.write(ofs, opt.k, sv_counter))
                    return false;
    }

    // I/O
    ofs.close();
    sc.write_CSV();     // OPTIONAL

    return true;
}


Traceback ExtendedCCDBG::DFS_Init_bidirectional(const UnitigColorMap<UnitigExtension> &ucm,
                                                const bool verbose,
                                                Setcover<std::unordered_set<unsigned int> > &sc){
    Traceback tb_bw;
    Traceback tb_fw;

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    if (verbose) cout << "Starting at " << ue->getID() << endl;

    BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();
    ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

    // ----------------
    // | go backward  |
    // ----------------
    if (!bw_neighbors.hasPredecessors()){
        tb_bw.recursive_return_status = true;
    }
    else{
        if (verbose) cout << "Setting " << ue->getID() << " to seen (bw)" << endl;
        ue->set_seen_bw();

        neighborsContainer descendingSortedNeighbors;
        sortNeighbors(ucm, bw_neighbors, descendingSortedNeighbors);

        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            if (tb_bw.recursive_return_status == false){    //restriction to only one successful path search
                for (auto &neighbor : bw_neighbors){
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;

                    if (ue_id == neighbor_id){
                        DFS_case_NEW(ucm, neighbor, tb_bw, sc, verbose);
                        break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                }
            }
            else {break;}
        }
    }

    // ---------------
    // | go forward  |
    // ---------------
    if (!fw_neighbors.hasSuccessors()){
        tb_fw.recursive_return_status = true;
    }
    else{
        if (verbose) cout << "Setting " << ue->getID() << " to seen (fw)" << endl;
        ue->set_seen_fw();

        neighborsContainer descendingSortedNeighbors;
        sortNeighbors(ucm, fw_neighbors, descendingSortedNeighbors);

        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            if (tb_fw.recursive_return_status == false){    //restriction to only one successful path search
                for (auto &neighbor : fw_neighbors){
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;

                    if (ue_id == neighbor_id){
                        DFS_case_NEW(ucm, neighbor, tb_fw, sc, verbose);
                        break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                }
            }
            else{break;}
        }
    }

    // ----------
    // | return |
    // ----------
    Traceback tb;

    if (tb_bw.recursive_return_status && tb_fw.recursive_return_status){                                                // both directions had successful DFS
        if (bw_neighbors.hasPredecessors() && fw_neighbors.hasSuccessors() && sc.hasMinContribution()){                 // started at internal node
            tb.recursive_return_status = true;
            tb.rearrange(tb_bw, tb_fw);
            sc.incorporate();
        }
        else if (bw_neighbors.hasPredecessors() && !fw_neighbors.hasSuccessors() && sc.hasMinContribution()){           // start node was an end node
            tb.recursive_return_status = true;
            tb.join(tb_bw);
            sc.incorporate();
        }
        else if (!bw_neighbors.hasPredecessors() && fw_neighbors.hasSuccessors() && sc.hasMinContribution()){           // start node was an end node
            tb.recursive_return_status = true;
            tb.join(tb_fw);
            sc.incorporate();
        }
        else if (!bw_neighbors.hasPredecessors() && !fw_neighbors.hasSuccessors()){                                     // start node is a singleton
            tb.recursive_return_status = true;
            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);
            sc.add(ue->getID());
            sc.incorporate();
        }
        else{
            cout << ue->getID() << " was rejected by the Setcover." << endl;
        }
    }
    else{
        cout << "At least one DFS traversal (bw/fw) was not successful." << endl;
    }

    if (verbose) cout << "Done with " << ue->getID() << endl;

    sc.clear();
    DFS_cleaner_seen_only();

    return tb;
}


Traceback ExtendedCCDBG::DFS_Visit_NEW(const UnitigColorMap<UnitigExtension> &ucm,
                                       const uint8_t src_direction,
                                       Setcover<> &sc,
                                       const bool verbose){
    Traceback tb;

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    uint8_t traversal_direction = (src_direction==GO_BACKWARD) ? GO_FORWARD : GO_BACKWARD;  // NOTE: this could be avoided if I'd pass traversal direction directly in DFS_case

    // --------------
    // |  backward  |
    // --------------
    if (traversal_direction==GO_BACKWARD){
        if (verbose) cout << "Setting " << ue->getID() << " seen (bw)." << endl;
        ue->set_seen_bw();

        // ------------------
        // |  if sink node  |
        // ------------------
        BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();

        if (!bw_neighbors.hasPredecessors()){
            if (verbose) cout << "Node " << ue->getID() << " is a sink. Start traceback." << endl;
            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);
            tb.recursive_return_status = true;
            sc.add(ue->getID());
            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------
        neighborsContainer descendingSortedNeighbors;
        sortNeighbors(ucm, bw_neighbors, descendingSortedNeighbors);

        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            if (tb.recursive_return_status == false){
                for (auto &neighbor : bw_neighbors){
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;
                    if (ue_id == neighbor_id){
                        DFS_case_NEW(ucm, neighbor, tb, sc, verbose);
                        break;      // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                }
            }
            else{break;}
        }
    }

    // -------------
    // |  forward  |
    // -------------
    else{   // traversal_direction==GO_FORWARD
        if (verbose) cout << "Setting " << ue->getID() << " seen (fw)." << endl;
        ue->set_seen_fw();

        // ------------------
        // |  if sink node  |
        // ------------------
        ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

        if (!fw_neighbors.hasSuccessors()){
            if (verbose) cout << "Node " << ue->getID() << " is a sink. Start traceback." << endl;
            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);
            tb.recursive_return_status = true;
            sc.add(ue->getID());
            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------
        neighborsContainer descendingSortedNeighbors;
        sortNeighbors(ucm, fw_neighbors, descendingSortedNeighbors);

        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            if (tb.recursive_return_status == false){
                for (auto &neighbor : fw_neighbors){
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;
                    if (ue_id == neighbor_id){
                        DFS_case_NEW(ucm, neighbor, tb, sc, verbose);
                        break;    // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                }
            }
            else{break;}
        }
    }

    return tb;
}


void ExtendedCCDBG::DFS_case_NEW(const UnitigColorMap<UnitigExtension> &ucm,
                                 const UnitigColorMap<UnitigExtension> &neighbor,
                                 Traceback &tb,
                                 Setcover<> &sc,
                                 const bool verbose){

    DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
    UnitigExtension* ucm_ue = ucm_da->getData(ucm);

    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);

    // ----------------
    // |   Case 1+3   |
    // ----------------
    if (whereFrom(neighbor, ucm)==GO_BACKWARD){
        /*  Case 3:                             Case 1:
        *   SRC ------->                OR              -------> SRC
        *           -------> SUC                PRE <--------
        */

        if ( neighbor_ue->is_undiscovered_fw() ){
            if (verbose) cout << "Traversal at " << ucm_ue->getID() << " will go to " << neighbor_ue->getID() << endl;
            sc.add(ucm_ue->getID());    // covers DFS_Init and DFS_Visit
            Traceback returned_tb = DFS_Visit_NEW(neighbor, GO_BACKWARD, sc, verbose);
            if (verbose) cout << "Jumped back to " << ucm_ue->getID() << endl;

            if (returned_tb.recursive_return_status){
                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);
            }
            else{
                sc.del(ucm_ue->getID());   // in case DFS jumps back without traceback, e.g. at loop
            }
        }
        else{
            // Unitig was seen before, happens e.g. at loops.
            if (verbose) cout << "Neighbor node " << neighbor_ue->getID() << " was seen before. Jump back without traceback." << endl;
        }
    }
    // ----------------
    // |   Case 2+4   |
    // ----------------
    else{   // whereFrom(neighbor, ucm)==GO_FORWARD
        /*  Case 4:                             Case 2:
        *   SRC ------->                OR              -------> SRC
        *           <------- SUC                PRE ------->
        */

        if ( neighbor_ue->is_undiscovered_bw() ){
            if (verbose) cout << "Traversal at " << ucm_ue->getID() << " will go to " << neighbor_ue->getID() << endl;
            sc.add(ucm_ue->getID());    // covers DFS_Init and DFS_Visit
            Traceback returned_tb = DFS_Visit_NEW(neighbor, GO_FORWARD, sc, verbose);
            if (verbose) cout << "Jumped back to " << ucm_ue->getID() << endl;

            if (returned_tb.recursive_return_status){
                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);
            }
            else{
                sc.del(ucm_ue->getID());   // in case DFS jumps back without traceback, e.g. at loop
            }
        }
        else{
            // Unitig was seen before, happens e.g. at loops.
            if (verbose) cout << "Neighbor node " << neighbor_ue->getID() << " was seen before. Jump back without traceback." << endl;
        }
    }

}


template <class TNeighborCDBG>
inline void ExtendedCCDBG::sortNeighbors(const UnitigColorMap<UnitigExtension> &ucm,
                                         const TNeighborCDBG &neighbors,
                                         neighborsContainer &container) const{
    //TODO: [Improvement] if cardinality of neighbor object is 1, check_common_colors() doesn't need to be called

    for (auto &neighbor : neighbors){
        DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
        UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
        unsigned neighbor_id = neighbor_ue->getID();
        unsigned cc_ = check_common_colors(ucm, neighbor);
        container.insert(std::pair<unsigned, unsigned>(cc_, neighbor_id));
    }
}


unsigned ExtendedCCDBG::check_common_colors(const UnitigColorMap <UnitigExtension> &ucm, const UnitigColorMap <UnitigExtension> &neighbor) const{
    const UnitigColors*      ucm_unitig_colors =      ucm.getData()->getUnitigColors(ucm);
    const UnitigColors* neighbor_unitig_colors = neighbor.getData()->getUnitigColors(ucm);

    size_t ucm_max_color_index = ucm_unitig_colors->colorMax(ucm);      //NOTE: min(colorMax(ucm),colorMax(neighbor)) for more efficiency?

    // compute #common colors in ucm and neighbor
    unsigned counter = 0;
    for (unsigned color_id=0; color_id<=ucm_max_color_index; ++color_id){
        if ((ucm_unitig_colors->size(ucm, color_id)!=0) && (neighbor_unitig_colors->size(neighbor, color_id!=0))){
            counter++;
        }
    }

    return counter;
}















