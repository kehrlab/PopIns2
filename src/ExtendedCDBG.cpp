#include "ExtendedCDBG.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>




// default constructor
ExtendedCDBG::ExtendedCDBG(int kmer_length, int minimizer_length) : CompactedDBG< UnitigExtension >(kmer_length, minimizer_length), init_status(false), isKmerCovInit(false), dfs_time(0), dfs_passed(false) {
    /* 1) IDs are not initiated at construction time (see init_ids())
     * 2) The UnionFind vector will be empty an construction time, will be resized at use.
     * 3) kmer_coverage vector will be empty an construction time, will be resized at use.
     */
}


void ExtendedCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &unitig : *this){
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


void ExtendedCDBG::print_unitig_info(){
    // -------------------------------------------------------------------------------------------
    // | ID | CC | Pre | Suc | DFS color | DFS predecessor| DFS discovery time | DFS finish time |
    // -------------------------------------------------------------------------------------------
    for (auto &unitig : *this){
        unsigned uid = unitig.getData()->getID();
        unsigned ucc = seqan::findSet(UF, uid);
        cout << "ID:" << uid << " | CC:" << ucc;

        BackwardCDBG<UnitigExtension, false> bw_dbg = unitig.getPredecessors();
        cout << " | Pre:";
        for (auto &predecessor : bw_dbg){
            size_t pre_id = predecessor.getData()->getID();
            cout << pre_id << ",";
        }

        ForwardCDBG<UnitigExtension, false> fw_dbg = unitig.getSuccessors();
        cout << " | Suc:";
        for (auto &successor : fw_dbg){
            size_t suc_id = successor.getData()->getID();
            cout << suc_id << ",";
        }

        cout << " | " << unitig.getData()->dfs_color << " | dfs-pre:" << unitig.getData()->dfs_ancestor << " | d.time:" << unitig.getData()->dfs_discovertime << " | f.time:" << unitig.getData()->dfs_finishtime;
        cout << endl;
    }
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
        size_t unitig_id = unitig.getData()->getID();

        BackwardCDBG<UnitigExtension, false> predecessors = unitig.getPredecessors();
        for (auto &it_pre : predecessors){
            unsigned pre_id = it_pre.getData()->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, pre_id));
        }

        ForwardCDBG<UnitigExtension, false> successors = unitig.getSuccessors();
        for (auto &it_suc : successors){
            unsigned suc_id = it_suc.getData()->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, suc_id));
        }
    }

    return true;
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


/*!
 * \fn      void ExtendedCDBG::dfs(UnitigMap<UnitigExtension> &um)
 * \brief   Local depth-first search (DFS) within a connected component, given a certain start node.
 */
void ExtendedCDBG::dfs(UnitigMap<UnitigExtension> &um){
    ForwardCDBG<UnitigExtension, false> successors = um.getSuccessors();
    for (auto &suc : successors)
        if (suc.getData()->dfs_color == 'w'){
            suc.getData()->dfs_ancestor = um.getData()->getID();
            dfs_visit(suc);
        }
    um.getData()->dfs_color = 'b';
    dfs_time += 1;
    um.getData()->dfs_finishtime = dfs_time;

    dfs_passed = true;
}


/*!
 * \fn      void ExtendedCDBG::dfs_visit()
 * \brief   Depth-first search neighbor traversal and visiting time updates
 */
void ExtendedCDBG::dfs_visit(UnitigMap<UnitigExtension> &um){
    dfs_time += 1;
    um.getData()->dfs_discovertime = dfs_time;
    um.getData()->dfs_color = 'g';

    // strand awareness
    // WARNING: I don't know whether this approach of strand awareness works yet.
    if (um.strand){     // if forward strand (+)
        ForwardCDBG<UnitigExtension, false> successors = um.getSuccessors();
        for (auto &suc : successors){
            if (suc.getData()->dfs_color == 'w'){
                suc.getData()->dfs_ancestor = um.getData()->getID();
                dfs_visit(suc);
            }
        }
    }
    else{               // if reverse complement (-)
        BackwardCDBG<UnitigExtension, false> predecessors = um.getPredecessors();
        for (auto &pre : predecessors){
            if (pre.getData()->dfs_color == 'w'){
                pre.getData()->dfs_ancestor = um.getData()->getID();
                dfs_visit(pre);
            }
        }
    }

    um.getData()->dfs_color = 'b';
    dfs_time += 1;
    um.getData()->dfs_finishtime = dfs_time;

    // if DFS hits a sink, return the current DFS path to source
    if (um.getData()->isSink()){
        std::vector<unsigned> path;
        path.push_back(um.getData()->getID());
        dfs_traceback(path, um);

        // [DEBUG] start
        for (unsigned id : path) cout << id << ",";
        cout << endl;
        // [DEBUG] end
    }
}


/*!
 * \fn          void ExtendedCDBG::dfs_traceback(vector<unsigned> &vec, UnitigMap<UnitigExtension> &um_sink)
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
void ExtendedCDBG::dfs_traceback(vector<unsigned> &vec, UnitigMap<UnitigExtension> &um_sink){

    unsigned dfs_anc = um_sink.getData()->dfs_ancestor;
    vec.push_back(dfs_anc);

    // if um_sink is a source node, then the dfs_ancestor is zero (default init) and we can stop traversing
    if (dfs_anc != 0){

        bool not_found_yet = true;

        BackwardCDBG<UnitigExtension, false> predecessors = um_sink.getPredecessors();
        for (auto &pre : predecessors)
            if (pre.getData()->getID() == dfs_anc){
                not_found_yet = false;
                dfs_traceback(vec, pre);
                break;
            }

        if (!not_found_yet) return;

        ForwardCDBG<UnitigExtension, false> successors = um_sink.getSuccessors();
        for (auto &suc : successors)
            if (suc.getData()->getID() == dfs_anc)
                dfs_traceback(vec, suc);

    }
}


/*!
 * \fn          bool ExtendedCDBG::annotate_kmer_coverage(const vector<string> &sample_fastx_names)
 * \brief       Reads the input files again to annotate a coverage per kmer
 * \return      true if successful
 */
// TODO: Runtime can potentially be improved with a specialized library/approach that is not loading in but just streaming kmers.
bool ExtendedCDBG::annotate_kmer_coverage(const vector<string> &sample_fastx_names){

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
                    UnitigMap<UnitigExtension> um = find(kmer);

                    // increment kmer count in unitig
                    if (um.size != 0)   // if kmer is in graph
                        if (um.getData()->kmer_coverage[um.dist] < 65535)   // if unsigned short saturation isn't reached yet
                            um.getData()->kmer_coverage[um.dist] += 1;

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
 * \fn      void ExtendedCDBG::init_kmer_cov()
 * \brief   Inits a tally to count the kmer coverage for each unitig.
 */
void ExtendedCDBG::init_kmer_cov(){
    for (auto &unitig : *this){
        size_t sz_v_kmer_cov = unitig.size - (this->getK() - 1);
        (unitig.getData()->kmer_coverage).resize(sz_v_kmer_cov);
    }
    isKmerCovInit = true;
}


/*!
 * \fn      void ExtendedCDBG::small_bubble_removal()
 * \brief   Remove bubbles within a certain max distance (in bases) from the start node.
 */
void ExtendedCDBG::small_bubble_removal(){
    size_t delta_k = getK()<<1;

    for (auto &unitig : *this){
        if (unitig.getData()->getID() == 257){      // TEST
            PathSet small_bubble_paths;
            if (!bfs_with_max_dist(unitig, small_bubble_paths, delta_k))
                cerr << "WARNING: an unexpected senario occured during the small bubble detection." << endl;
            clear_path_search_attributes();

            cout << "----------" << endl;       // TEST
            for (auto &set : small_bubble_paths){
                cout << "[";
                for (auto &um : set){
                    cout << um.getData()->getID() << ", ";
                }
                cout << "]" << endl;
            }
            cout << "----------" << endl;       // TEST
        }
    }
}


/*!
 * \fn      void ExtendedCDBG::bfs_with_max_dist()
 * \brief   Run a BFS until the max distance in bases from the startnode (um) is exceeded.
 * \return  true if successful
 */
inline bool ExtendedCDBG::bfs_with_max_dist(UnitigMap<UnitigExtension> &um, PathSet &pathset, const size_t max_dist){

    um.getData()->dfs_color = 's';  // 's' is special color of start node

    std::queue<UnitigMap<UnitigExtension>> q;

    q.push(um);

    while (!q.empty()){

        UnitigMap<UnitigExtension> um_ = q.front();
        q.pop();

        for (auto &pre : um_.getPredecessors()){

            // find undiscovered node in BFS < max_dist, push node in queue
            if (pre.getData()->dfs_color == 'w' && um_.getData()->dfs_discovertime + pre.size < max_dist){
                pre.getData()->dfs_color = 'g';
                pre.getData()->dfs_discovertime = um_.getData()->dfs_discovertime + pre.size;
                pre.getData()->dfs_ancestor = um_.getData()->getID();
                q.push(pre);
            }

            // find undiscovered node in BFS >= max_dist, break traversal
            else if (pre.getData()->dfs_color == 'w' && um_.getData()->dfs_discovertime + pre.size >= max_dist){
                pre.getData()->dfs_color = 'g';
            }

            // if already discovered node is >= max_dist (bases) away from source, we found a 'small bubble'
            else if (pre.getData()->dfs_color == 'g' && um_.getData()->dfs_discovertime + pre.size >= max_dist && pre.getData()->dfs_color != 's'){
                UnitigPath up;
                if (get_reverse_bfs_paths(pre, up)){
                    pathset.push_back(up);
                }
                else{
                    cerr << "BFS traceback returned an error." << endl;
                    return 0;
                }
            }

            // BFS returned to start node
            else if (pre.getData()->dfs_color == 's'){
                continue;
            }

            else{
                cerr << "WARNING: BFS (pre) ran into an undefined case from " << um_.getData()->getID() << "(um) to " << pre.getData()->getID() << "(pre)." << endl;
                return 0;
            }
        }

        for (auto &suc : um_.getSuccessors()){

            // find undiscovered nodes in BFS < max_dist, push node in queue
            if (suc.getData()->dfs_color == 'w' && um_.getData()->dfs_discovertime + suc.size < max_dist){
                suc.getData()->dfs_color = 'g';
                suc.getData()->dfs_discovertime = um_.getData()->dfs_discovertime + suc.size;
                suc.getData()->dfs_ancestor = um_.getData()->getID();
                q.push(suc);
            }

            // find undiscovered node in BFS >= max_dist, break traversal
            else if (suc.getData()->dfs_color == 'w' && um_.getData()->dfs_discovertime + suc.size >= max_dist){
                suc.getData()->dfs_color = 'g';
            }

            // if already discovered node is >= max_dist (bases) away from source, we found a 'small bubble'
            else if (suc.getData()->dfs_color == 'g' && um_.getData()->dfs_discovertime + suc.size >= max_dist && suc.getData()->dfs_color != 's'){
                UnitigPath up;
                if (get_reverse_bfs_paths(suc, up)){
                    pathset.push_back(up);
                }
                else{
                    cerr << "BFS traceback returned an error." << endl;
                    return 0;
                }

            }

            // BFS returned to start node
            else if (suc.getData()->dfs_color == 's'){
                continue;
            }

            else{
                cerr << "WARNING: BFS (suc) ran into an undefined case from " << um_.getData()->getID() << "(um) to " << suc.getData()->getID() << "(suc)." << endl;
                return 0;
            }

        }
    }

    return 1;
}


/*!
 * \fn      bool ExtendedCDBG::get_reverse_bfs_paths(UnitigMap<UnitigExtension> &um)
 * \brief   This function returns all the bubble paths from sink to source.
 * \return  true if successful
 */
bool ExtendedCDBG::get_reverse_bfs_paths(UnitigMap<UnitigExtension> &um, UnitigPath &up){

    bool ret_ = false;
    bool barrier_ = false;

    for (auto &pre : um.getPredecessors()){
        // path node
        if (pre.getData()->dfs_ancestor != 0){
            up.push_back(pre);
            ret_ = get_reverse_bfs_paths(pre, up);
        }

        // found node outside the max_dist area of BFS
        else if (pre.getData()->dfs_ancestor == 0 && pre.getData()->dfs_color == 'w' && pre.getData()->dfs_color != 's'){
            continue;
        }

        // former traceback source found
        else if (pre.getData()->dfs_ancestor == 0 && pre.getData()->dfs_color == 'g' && pre.getData()->dfs_color != 's'){
            continue;
        }

        // former BFS source found
        else if (pre.getData()->dfs_color == 's'){
            barrier_ = true;
            ret_ = true;
            break;
        }

        else{
            cerr << "ERROR: Small bubble-popping traceback ran into an undefined case from " << um.getData()->getID() << "(um) to " << pre.getData()->getID() << "(pre) while tracing predecessors." << endl;
            return ret_;
        }
    }

    if (barrier_) return ret_;

    for (auto &suc : um.getSuccessors()){
        // path node
        if (suc.getData()->dfs_ancestor != 0){
            up.push_back(suc);
            ret_ = get_reverse_bfs_paths(suc, up);
        }

        // found node outside the max_dist area of BFS
        else if (suc.getData()->dfs_ancestor == 0 && suc.getData()->dfs_color == 'w' && suc.getData()->dfs_color != 's'){
            continue;
        }

        // former traceback source found
        else if (suc.getData()->dfs_ancestor == 0 && suc.getData()->dfs_color == 'g' && suc.getData()->dfs_color != 's'){
            continue;
        }

        // former BFS source found
        else if (suc.getData()->dfs_color == 's'){
            ret_ = true;
            break;
        }

        else{
            cerr << "ERROR: Small bubble-popping traceback ran into an undefined case from " << um.getData()->getID() << "(um) to " << suc.getData()->getID() << "(suc) while tracing successors." << endl;
            return ret_;
        }
    }

    return ret_;
}


/*!
 * \fn      void ExtendedCDBG::clear_path_search_attributes()
 * \brief   This functions resets all xfs (dfs/bfs) attributes to its default.
 */
inline void ExtendedCDBG::clear_path_search_attributes(){
    for (auto &unitig : *this){
        unitig.getData()->dfs_color = 'w';
        unitig.getData()->dfs_ancestor = 0;
        unitig.getData()->dfs_discovertime = 0;
        unitig.getData()->dfs_finishtime = 0;
    }
}



































