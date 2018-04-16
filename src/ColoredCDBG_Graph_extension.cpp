#include "ColoredCDBG_Graph_extension.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>




// default constructor
ExtendedCCDBG::ExtendedCCDBG(int kmer_length, int minimizer_length) : ColoredCDBG<UnitigExtension> (kmer_length, minimizer_length), init_status(false), isKmerCovInit(false), dfs_time(0), dfs_passed(false) {
    /* 1) IDs are not initiated at construction time (see init_ids())
     * 2) The UnionFind vector will be empty an construction time, will be resized at use.
     * 3) kmer_coverage vector will be empty an construction time, will be resized at use.
     */
}

/*
inline UnitigExtension* ExtendedCCDBG::getDataDirectly(UnitigColorMap &ucm){
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* data = da->getData(ucm);
    return data;
}
*/

void ExtendedCCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);      // ue is a POINTER to a UnitigExtension
        ue->setID(i);
        ++i;
    }
    init_status = true;
}


void ExtendedCCDBG::print_ids(){
    if (is_init() == true){
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


void ExtendedCCDBG::print_unitig_info(){
    // -------------------------------------------------------------------------------------------
    // | ID | CC | Pre | Suc | DFS color | DFS predecessor| DFS discovery time | DFS finish time |
    // -------------------------------------------------------------------------------------------
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        unsigned uid = ue->getID();
        unsigned ucc = seqan::findSet(UF, uid);
        cout << "ID:" << uid << " | CC:" << ucc;

        cout << " | Pre:";
        for (auto &predecessor : unitig.getPredecessors()){
            DataAccessor<UnitigExtension>* da = predecessor.getData();
            UnitigExtension* ue = da->getData(predecessor);
            size_t pre_id = ue->getID();
            cout << pre_id << ",";
        }

        cout << " | Suc:";
        for (auto &successor : unitig.getSuccessors()){
            DataAccessor<UnitigExtension>* da = successor.getData();
            UnitigExtension* ue = da->getData(successor);
            size_t suc_id = ue->getID();
            cout << suc_id << ",";
        }

        cout << " | " << ue->dfs_color << " | dfs-pre:" << ue->dfs_ancestor << " | d.time:" << ue->dfs_discovertime << " | f.time:" << ue->dfs_finishtime;
        cout << endl;
    }
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


/*!
 * \fn      void ExtendedCCDBG::dfs(UnitigColorMap<UnitigExtension> &um)
 * \brief   Local depth-first search (DFS) within a connected component, given a certain start node.
 */
void ExtendedCCDBG::dfs(const UnitigColorMap<UnitigExtension> &um){
    DataAccessor<UnitigExtension>* da = um.getData();
    UnitigExtension* um_ue = da->getData(um);

    for (auto &suc : um.getSuccessors()){
        DataAccessor<UnitigExtension>* da = suc.getData();
        UnitigExtension* ue = da->getData(suc);
        if (ue->dfs_color == 'w'){
            ue->dfs_ancestor = um_ue->getID();
            dfs_visit(suc);
        }
    }
    um_ue->dfs_color = 'b';
    dfs_time += 1;
    um_ue->dfs_finishtime = dfs_time;

    dfs_passed = true;
}


/*!
 * \fn      void ExtendedCCDBG::dfs_visit(UnitigColorMap<UnitigExtension> &um)
 * \brief   Depth-first search neighbor traversal and visiting time updates
 */
void ExtendedCCDBG::dfs_visit(const UnitigColorMap<UnitigExtension> &um){
    dfs_time += 1;
    DataAccessor<UnitigExtension>* um_da = um.getData();
    UnitigExtension* um_ue = um_da->getData(um);
    um_ue->dfs_discovertime = dfs_time;
    um_ue->dfs_color = 'g';

    // strand awareness
    // WARNING: I don't know whether this approach of strand awareness works yet.
    if (um.strand){     // if forward strand (+)
        for (auto &suc : um.getSuccessors()){
            DataAccessor<UnitigExtension>* suc_da = suc.getData();
            UnitigExtension* suc_ue = suc_da->getData(suc);
            if (suc_ue->dfs_color == 'w'){
                suc_ue->dfs_ancestor = um_ue->getID();
                dfs_visit(suc);
            }
        }
    }
    else{               // if reverse complement (-)
        for (auto &pre : um.getPredecessors()){
            DataAccessor<UnitigExtension>* pre_da = pre.getData();
            UnitigExtension* pre_ue = pre_da->getData(pre);
            if (pre_ue->dfs_color == 'w'){
                pre_ue->dfs_ancestor = um_ue->getID();
                dfs_visit(pre);
            }
        }
    }

    um_ue->dfs_color = 'b';
    dfs_time += 1;
    um_ue->dfs_finishtime = dfs_time;

    // if DFS hits a sink, return the current DFS path to source
    if (um_ue->isSink()){
        std::vector<unsigned> path;
        path.push_back(um_ue->getID());
        dfs_traceback(path, um);

        // [DEBUG] start
        for (unsigned id : path) cout << id << ",";
        cout << endl;
        // [DEBUG] end
    }
}


/*!
 * \fn          void ExtendedCCDBG::dfs_traceback(vector<unsigned> &vec, const UnitigMap<UnitigExtension> &um_sink)
 * \brief       This function walks back on the DFS trace of the nodes, from a given sink to its source.
 * \details     It is important for the dfs_traceback to know in which orientation the unitigs are. This was either
 *              predecessors or ancestors have to be choosen for further path traversal. E.g.:
 *                 -----> 1
 *                    -----> 2
 *                        <----- rc(3)
 *                           <----- rc(4)
 *                               -----> 5
 *              Trace: 5, rc(4), rc(3), 2, 1
 * \return      a list containing the path from sink to source, i.e. in reverse directional order
 */
/* TODO: I am checking too many neighbors here, i.e. I do up to 8 comparisons while I only had to
 * do up to 4, because I couldn't figure out how to trace RC properly. But it works.
 */
void ExtendedCCDBG::dfs_traceback(vector<unsigned> &vec, const UnitigColorMap<UnitigExtension> &um_sink){
    DataAccessor<UnitigExtension>* da = um_sink.getData();
    UnitigExtension* ue = da->getData(um_sink);

    unsigned dfs_anc = ue->dfs_ancestor;
    vec.push_back(dfs_anc);

    // if um_sink is a source node, then the dfs_ancestor is zero (default init) and we can stop traversing
    if (dfs_anc != 0){

        bool not_found_yet = true;

        for (auto &pre : um_sink.getPredecessors()){
            DataAccessor<UnitigExtension>* pre_da = pre.getData();
            UnitigExtension* pre_ue = pre_da->getData(pre);
            if (pre_ue->getID() == dfs_anc){
                not_found_yet = false;
                dfs_traceback(vec, pre);
                break;
            }
        }

        if (!not_found_yet) return;

        for (auto &suc : um_sink.getSuccessors()){
            DataAccessor<UnitigExtension>* suc_da = suc.getData();
            UnitigExtension* suc_ue = suc_da->getData(suc);
            if (suc_ue->getID() == dfs_anc)
                dfs_traceback(vec, suc);
        }

    }
}


/*!
 * \fn          bool ExtendedCCDBG::annotate_kmer_coverage(const vector<string> &sample_fastx_names)
 * \brief       Reads the input files again to annotate a coverage per kmer
 * \return      true if successful
 */
// TODO: Runtime can potentially be improved with a specialized library/approach that is not loading in but just streaming kmers.
bool ExtendedCCDBG::annotate_kmer_coverage(const vector<string> &sample_fastx_names){

    if (!isKmerCovInit){
        cerr << "ERROR: kmer_coverage has not been initialized." << endl;
        return 0;
    }

    // reading one file at a time
    std::vector<string>::const_iterator itFile = sample_fastx_names.cbegin();
    while (itFile != sample_fastx_names.cend()){

        seqan::CharString seqFileName = seqan::toCString(*itFile);
        seqan::SeqFileIn seqFileIn;

        if (!seqan::open(seqFileIn, seqan::toCString(seqFileName))){
            // if file could not be found, retry with "./" directory prefix
            seqan::CharString pref = "./";
            seqan::append(pref, seqFileName);
            if (!seqan::open(seqFileIn, seqan::toCString(pref))){
                std::cerr << "ERROR: Could not open the file.\n";
                return 0;
            }
        }

        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::DnaString> seqs;
        try{
            seqan::readRecords(ids, seqs, seqFileIn);
            // iterate reads
            typedef seqan::Iterator<seqan::StringSet<seqan::DnaString> >::Type TStringSetIterator;
            for (TStringSetIterator itSeq = seqan::begin(seqs); itSeq != seqan::end(seqs); ++itSeq){
                // iterate kmers
                size_t endpos = seqan::endPosition(*itSeq) - getK();
                for (size_t pos = 0; pos <= endpos; ++pos){
                    seqan::String<char, seqan::CStyle> infix = seqan::infixWithLength(*itSeq, pos, getK());
                    Kmer kmer(infix);
                    UnitigColorMap<UnitigExtension> um = find(kmer);
                    DataAccessor<UnitigExtension>* da = um.getData();
                    UnitigExtension* ue = da->getData(um);

                    // increment kmer count in unitig
                    if (um.size != 0)   // if kmer is in graph
                        if (ue->kmer_coverage[um.dist] < 65535)   // if unsigned short saturation isn't reached yet
                            ue->kmer_coverage[um.dist] += 1;

                    // NOTE: The function finds all kmers. However, some of the following UnitigMaps have length 0.
                    //       The only explanation to me is that the corresponding unitigs got filtered out. Since
                    //       filter() erases unitigs and does NOT leave a mask I cannot tell 100%.
                    //       There are all coverage vectors filled without any 'zero gaps'. This indicates that
                    //       the coverage was successfully recovered. The minimum kmer coverage is 2, where forward
                    //       and RC count both for the same kmer.
                }
            }
        }
        catch (exception const & e){
            std::cout << "ERROR: " << e.what() << std::endl;
            return 0;
        }

        seqan::close(seqFileIn);
        itFile++;
    }
    return 1;
}


/*!
 * \fn      void ExtendedCCDBG::init_kmer_cov()
 * \brief   Inits a tally to count the kmer coverage for each unitig.
 */
void ExtendedCCDBG::init_kmer_cov(){
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);

        size_t sz_v_kmer_cov = unitig.size - (this->getK() - 1);
        (ue->kmer_coverage).resize(sz_v_kmer_cov);
    }
    isKmerCovInit = true;
}


/*!
 * \fn      void ExtendedCCDBG::small_bubble_removal()
 * \brief   Remove bubbles within a certain max distance (in bases) from the start node.
 */
void ExtendedCCDBG::small_bubble_removal(){
    size_t delta_k = getK()<<1;

    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);

        //if (unitig.getData()->getID() == 257){      // TEST
            PathSet small_bubble_paths;
            if (!bfs_with_max_dist(unitig, small_bubble_paths, delta_k))
                cerr << "WARNING: an unexpected senario occured during the small bubble detection." << endl;
            clear_path_search_attributes();

            cout << "----" << ue->getID() << "----" << endl;       // TEST
            cout << "[";
            for (auto &set : small_bubble_paths){
                cout << "[";
                for (auto &um : set){
                    DataAccessor<UnitigExtension>* um_da = um.getData();
                    UnitigExtension* um_ue = um_da->getData(um);
                    cout << um_ue->getID() << ", ";
                }
                cout << "]";
            }
            cout << "]" << endl;
            cout << "----------" << endl << endl;       // TEST
        //}
    }
}


/*!
 * \fn      void ExtendedCCDBG::bfs_with_max_dist(const UnitigColorMap<UnitigExtension> &um, PathSet &pathset, const size_t max_dist)
 * \brief   Run a BFS until the max distance in bases from the startnode (um) is exceeded.
 * \return  true if successful
 */
inline bool ExtendedCCDBG::bfs_with_max_dist(const UnitigColorMap<UnitigExtension> &um, PathSet &pathset, const size_t max_dist){
    DataAccessor<UnitigExtension>* da = um.getData();
    UnitigExtension* ue = da->getData(um);

    ue->dfs_color = 's';  // 's' is special color of start node

    std::queue<UnitigColorMap<UnitigExtension>> q;

    q.push(um);

    while (!q.empty()){

        UnitigColorMap<UnitigExtension> um_ = q.front();
        q.pop();

        DataAccessor<UnitigExtension>* da_ = um_.getData();
        UnitigExtension* ue_ = da_->getData(um_);

        for (auto &pre : um_.getPredecessors()){

            DataAccessor<UnitigExtension>* pre_da = pre.getData();
            UnitigExtension* pre_ue = pre_da->getData(pre);

            // find undiscovered node in BFS < max_dist, push node in queue
            if (pre_ue->dfs_color == 'w' && ue_->dfs_discovertime + pre.size < max_dist){
                pre_ue->dfs_color = 'g';
                pre_ue->dfs_discovertime = ue_->dfs_discovertime + pre.size;
                pre_ue->dfs_ancestor = ue_->getID();
                q.push(pre);
            }

            // find undiscovered node in BFS >= max_dist, break traversal
            else if (pre_ue->dfs_color == 'w' && ue_->dfs_discovertime + pre.size >= max_dist){
                pre_ue->dfs_color = 'g';
            }

            // if already discovered node is >= max_dist (bases) away from source, we possibly found a small bubble
            else if (pre_ue->dfs_color == 'g' && ue_->dfs_discovertime + pre.size >= max_dist && pre_ue->dfs_color != 's'){
                // small bubble detected
                if (ue_->getID() != pre_ue->getID()){
                    UnitigPath up;
                    if (get_reverse_bfs_paths(pre, up)){
                        if (!up.empty()){
                            pathset.push_back(up);
                        }
                        else{
                            cerr << "CATCH CASE: [7]: Selfloop detected." << endl;
                            break;
                        }
                    }
                    else{
                        cerr << "WARNING: BFS traceback returned an error." << endl;
                        return 0;
                    }
                }
                else if (ue_->getID() == pre_ue->getID()){
                    cerr << "CATCH CASE: [5]: Selfloop detected." << endl;
                    break;      // segmentation fault (core dump) if I put it only 'continue'
                }
                else{
                    cerr << "WARNING: [1]: BFS (pre) ran into an undefined case from " << ue_->getID() << "(um) to " << pre_ue->getID() << "(pre)." << endl;
                    return 0;
                }
            }

            // BFS returned to start node
            else if (pre_ue->dfs_color == 's'){
                continue;
            }

            else{
                cerr << "WARNING: [2]: BFS (pre) ran into an undefined case from " << ue_->getID() << "(um) to " << pre_ue->getID() << "(pre)." << endl;
                return 0;
            }
        }

        for (auto &suc : um_.getSuccessors()){

            DataAccessor<UnitigExtension>* suc_da = suc.getData();
            UnitigExtension* suc_ue = suc_da->getData(suc);

            // find undiscovered nodes in BFS < max_dist, push node in queue
            if (suc_ue->dfs_color == 'w' && ue_->dfs_discovertime + suc.size < max_dist){
                suc_ue->dfs_color = 'g';
                suc_ue->dfs_discovertime = ue_->dfs_discovertime + suc.size;
                suc_ue->dfs_ancestor = ue_->getID();
                q.push(suc);
            }

            // find undiscovered node in BFS >= max_dist, break traversal
            else if (suc_ue->dfs_color == 'w' && ue_->dfs_discovertime + suc.size >= max_dist){
                suc_ue->dfs_color = 'g';
            }

            // if already discovered node is >= max_dist (bases) away from source, we possibly found a 'small bubble'
            else if (suc_ue->dfs_color == 'g' && ue_->dfs_discovertime + suc.size >= max_dist && suc_ue->dfs_color != 's'){
                // small bubble detected
                if (ue_->getID() != suc_ue->getID()){
                    UnitigPath up;
                    if (get_reverse_bfs_paths(suc, up)){
                        if (!up.empty()){
                            pathset.push_back(up);
                        }
                        else{
                            cerr << "CATCH CASE: [8]: Selfloop detected." << endl;
                            break;
                        }
                    }
                    else{
                        cerr << "WARNING: BFS traceback returned an error." << endl;
                        return 0;
                    }
                }
                else if (ue_->getID() == suc_ue->getID()){
                    cerr << "CATCH CASE: [6]: Selfloop detected." << endl;
                    break;      // segmentation fault (core dump) if I put it only 'continue'
                }
                else{
                    cerr << "WARNING: [3]: BFS (suc) ran into an undefined case from " << ue_->getID() << "(um) to " << suc_ue->getID() << "(suc)." << endl;
                    return 0;
                }
            }

            // BFS returned to start node
            else if (suc_ue->dfs_color == 's'){
                continue;
            }

            // small ring/circle/circuit detected, NOTE this code only needs to be in the successors loop since the predecessors need to mark the partner unitig first
            else if (suc_ue->dfs_color == 'g' && ue_->dfs_discovertime + suc.size < max_dist && suc_ue->dfs_color != 's'){
                // What TODO with this case?
                continue;
            }

            else{
                cerr << "WARNING: [4]: BFS (suc) ran into an undefined case from " << ue_->getID() << "(um) to " << suc_ue->getID() << "(suc)." << endl;
                return 0;
            }

        }
    }

    return 1;
}


/*!
 * \fn      bool ExtendedCCDBG::get_reverse_bfs_paths(const UnitigColorMap<UnitigExtension> &um, UnitigPath &up)
 * \brief   This function returns all the small bubble paths from sink to source.
 * \return  true if successful
 */
bool ExtendedCCDBG::get_reverse_bfs_paths(const UnitigColorMap<UnitigExtension> &um, UnitigPath &up){
    DataAccessor<UnitigExtension>* da = um.getData();
    UnitigExtension* ue = da->getData(um);

    bool ret_ = false;
    bool barrier_ = false;

    for (auto &pre : um.getPredecessors()){
        DataAccessor<UnitigExtension>* pre_da = pre.getData();
        UnitigExtension* pre_ue = pre_da->getData(pre);

        // path node
        if (pre_ue->dfs_ancestor != 0){
            up.push_back(pre);
            ret_ = get_reverse_bfs_paths(pre, up);
        }

        // found node outside the max_dist area of BFS
        else if (pre_ue->dfs_ancestor == 0 && pre_ue->dfs_color == 'w' && pre_ue->dfs_color != 's'){
            continue;
        }

        // former traceback source found
        else if (pre_ue->dfs_ancestor == 0 && pre_ue->dfs_color == 'g' && pre_ue->dfs_color != 's'){
            continue;
        }

        // former BFS source found
        else if (pre_ue->dfs_color == 's'){
            barrier_ = true;
            ret_ = true;
            break;
        }

        else{
            cerr << "ERROR: [9]: Small bubble-popping traceback ran into an undefined case from " << ue->getID() << "(um) to " << pre_ue->getID() << "(pre) while tracing predecessors." << endl;
            return ret_;
        }
    }

    if (barrier_) return ret_;

    for (auto &suc : um.getSuccessors()){
        DataAccessor<UnitigExtension>* suc_da = suc.getData();
        UnitigExtension* suc_ue = suc_da->getData(suc);
        // path node
        if (suc_ue->dfs_ancestor != 0){
            up.push_back(suc);
            ret_ = get_reverse_bfs_paths(suc, up);
        }

        // found node outside the max_dist area of BFS
        else if (suc_ue->dfs_ancestor == 0 && suc_ue->dfs_color == 'w' && suc_ue->dfs_color != 's'){
            continue;
        }

        // former traceback source found
        else if (suc_ue->dfs_ancestor == 0 && suc_ue->dfs_color == 'g' && suc_ue->dfs_color != 's'){
            continue;
        }

        // former BFS source found
        else if (suc_ue->dfs_color == 's'){
            ret_ = true;
            break;
        }

        else{
            cerr << "ERROR: [10]: Small bubble-popping traceback ran into an undefined case from " << ue->getID() << "(um) to " << suc_ue->getID() << "(suc) while tracing successors." << endl;
            return ret_;
        }
    }

    return ret_;
}


/*!
 * \fn      void ExtendedCCDBG::clear_path_search_attributes()
 * \brief   This functions resets all xfs (dfs/bfs) attributes to its default.
 */
inline void ExtendedCCDBG::clear_path_search_attributes(){
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);

        ue->dfs_color = 'w';
        ue->dfs_ancestor = 0;
        ue->dfs_discovertime = 0;
        ue->dfs_finishtime = 0;
    }
}



































