#include "ExtendedCDBG.h"
#include <unordered_set>
#include <unordered_map>



// default constructor
ExtendedCDBG::ExtendedCDBG(int kmer_length, int minimizer_length) : CompactedDBG< UnitigExtension >(kmer_length, minimizer_length), init_status(false), dfs_time(0), dfs_passed(false) {
    /* 1) IDs are not initiated at construction time (see init_ids())
     * 2) The UnionFind vector will be empty an construction time, will be resized at use.
     */
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


void ExtendedCDBG::print_unitig_info(){
    // --------------------------------------------------------------------------
    // | ID | CC | Pre | Suc | DFS color | DFS discovery time | DFS finish time |
    // --------------------------------------------------------------------------
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

        cout << " | " << unitig.getData()->dfs_color << " | d.time:" << unitig.getData()->dfs_discovertime << " | f.time:" << unitig.getData()->dfs_finishtime;
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
 * \fn      void ExtendedCDBG::dfs()
 * \brief   Depth-first search
 */
void ExtendedCDBG::dfs(){
    for (auto &unitig : *this)
        if (unitig.getData()->dfs_color == 'w')
            dfs_visit(unitig);

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

    ForwardCDBG<UnitigExtension, false> successors = um.getSuccessors();
    for (auto &suc : successors){
        if (suc.getData()->dfs_color == 'w'){
            suc.getData()->dfs_ancestor = um.getData()->getID();
            dfs_visit(suc);
        }
    }

    um.getData()->dfs_color = 'b';
    dfs_time += 1;
    um.getData()->dfs_finishtime = dfs_time;
}


/*!
 * \fn          bool ExtendedCDBG::single_source_longest_paths(UnitigMap<UnitigExtension> &um)
 * \brief       Traverse the graph and find all longest paths from a given start node.
 * \return      true if successful
 */
bool ExtendedCDBG::single_source_longest_paths(UnitigMap<UnitigExtension> &um){
    if (!dfs_passed)
        return 0;

    unsigned uid = um.getData()->getID();
    unsigned ucc = seqan::findSet(UF, uid);

    /*  FIXME: it is an incredible overhead to check ALL nodes here again whether they are in the same component.
     *  I need a data structure within ExtendedCDBG so that I only need to work on unitigs of the same CC.
     */
    // find all sinks in the connected component ("unitig" is a potential sink in this loop)
    for (auto &unitig : *this){
        // only continue within same connected component
        if (ucc != seqan::findSet(UF, unitig.getData()->getID()))
            continue;

        // check if sink
        if (unitig.getData()->dfs_finishtime - unitig.getData()->dfs_discovertime == 1){
            std::vector<unsigned> path;
            path.push_back(unitig.getData()->getID());
            traceback(path, unitig);

            // [DEBUG] start
            for (unsigned id : path) cout << id << ",";
            cout << endl;
            // [DEBUG] end
        }
    }



    return 1;
}


/*!
 * \fn          void ExtendedCDBG::traceback(vector<unsigned> &vec, UnitigMap<UnitigExtension> &um_sink)
 * \brief       This function walks back on the DFS trace of the nodes, from a given sink to its source.
 * \return      a list containing the path from sink to source, i.e. in reverse directional order
 */
void ExtendedCDBG::traceback(vector<unsigned> &vec, UnitigMap<UnitigExtension> &um_sink){
    unsigned dfs_anc = um_sink.getData()->dfs_ancestor;
    vec.push_back(dfs_anc);

    // if um_sink is a source node, then the dfs_ancestor is zero (default init) and we can stop traversing
    if (dfs_anc != 0){
        /* The only shorter way than looking for the ID in the whole graph again is by using the fact that the node with
        * the ancestor ID must be within the four predecessors of the unitig.
        */
        BackwardCDBG<UnitigExtension, false> predecessors = um_sink.getPredecessors();
        for (auto &pre : predecessors)
            if (pre.getData()->getID() == dfs_anc)
                traceback(vec, pre);
    }
}














