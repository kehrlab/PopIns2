#include "CompactedDBG_Graph_extension.h"



void ExtendedCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &um : *this){
        DataExtension* de = um.getData();      // de is a POINTER to a UnitigExtension
        de->setID(i);
        ++i;
    }
    id_init_status = true;
}


inline uint8_t ExtendedCDBG::getDegree(const UnitigMap<DataExtension> &um) const{
    uint8_t ret = 0;
    for (auto p : um.getPredecessors()) if (!p.isEmpty) ret++;
    for (auto s : um.getSuccessors())   if (!s.isEmpty) ret++;
    return ret;
}


/*!
 * \fn      const uint8_t ExtendedCDBG::whereToGo(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src) const
 * \brief   This function tests the predecessors P of a unitig u. If the searched unitig src is in P, then
 *          whereToGo() returns GO_FORWARD, denoting the traversal has to continue in the successors of u. If src is
 *          not in P, then whereToGo() returns GO_BACKWARD, denoting the traversal has to continue in the predecessors
 *          of u.
 * \details NOTE: This function is a lowest-level indication where to go, independent of the unitig's orientation. If
 *          the orientation of a unitig is rev-comp, then the result might be GO_BACKWARD while we still consider it
 *          a forward motion with respect to the traversal.
 * \return  DIRECTION
 */
inline uint8_t ExtendedCDBG::whereToGo(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src) const{
    uint8_t ret = GO_BACKWARD;
    for (auto &predecessor : um.getPredecessors())
         if (predecessor == src)    /* TODO: if small-loop it will always evaluate to true */
             ret = GO_FORWARD; 
    return ret;
}


/*!
 * \fn      bool ExtendedCDBG::DFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction)
 * \brief   This function initiates the recursion of the directed DFS.
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::DFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction){

    bool ret = 0;
    unsigned int dist = getK()-1;

    // mark current unitig
    DataExtension* de = um.getData();
    de->set_visited();

    um.strand==1? cout << "[[" << de->getID() << "+" << "]]" << endl : cout << "[[" << de->getID() << "-" << "]]" << endl;

    if (direction==GO_BACKWARD){
        for (auto &predecessor : um.getPredecessors()){
            DataExtension* de_pre = predecessor.getData();
            if (predecessor.size < DFS_MAX_DIST){
                cout << "\t[  ] I am at " << de->getID() << " and will go backward to " << de_pre->getID() << endl;
                /*  According to Stack Overflow [https://stackoverflow.com/questions/9021049/operator-for-a-boolean-in-c],
                *  the rhs of |= is always evaluated. That's wanted here since I don't want to break the recursion
                *  in new undefined cases.
                */
                ret |= DFS_Direction_Recursion(predecessor, um, dist);
                cout << "\tI jumped back to ID " << de->getID() << endl;
            }
        }
    }
    else if (direction==GO_FORWARD){
        for (auto &successor : um.getSuccessors()){
            DataExtension* de_suc = successor.getData();
            if (successor.size < DFS_MAX_DIST){
                cout << "\t[  ] I am at " << de->getID() << " and will go forward to " << de_suc->getID() << endl;
                /* see predecessor for equivalent explanation */
                ret |= DFS_Direction_Recursion(successor, um, dist);
                cout << "\tI jumped back to ID " << de->getID() << endl;
            }
        }
    }
    else{
        cerr << "ERROR: " << (direction) << " is not a valid direction" << endl;
        return 1;
    }

    return ret;
}


/*!
 * \fn      bool ExtendedCDBG::DFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist)
 * \brief   This function defines the recursion of the directed DFS.
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::DFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist){

    bool ret = 0;
    dist += um.size-(getK()-1);

    // unitig is of interest
    DataExtension* de = um.getData();

    if (dist>=DFS_MAX_DIST && de->is_not_visited()){
        cout << "\t[ 1] case exceeds max distance here at ID " << de->getID() << endl;
        de->set_visited();
        return ret;
    }

    um.strand==1? cout << "\t[[" << de->getID() << "+" << "]]" << endl : cout << "\t[[" << de->getID() << "-" << "]]" << endl;

    if (whereToGo(um, src)==GO_BACKWARD){
        for (auto &predecessor : um.getPredecessors()){
            DataExtension* de_pre = predecessor.getData();

            // NOTE: case order matters atm!

            if (predecessor==src){
                cout << "\t[ 2] small-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (predecessor==um){
                cout << "\t[ 3] self-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (de_pre->is_visited()){
                cout << "\t[ 4] small bubble detected ending at ID " << de_pre->getID() << endl;
                de->set_visited();
                return Traceback_Init(um, predecessor);
            }
            /* NOTE: What is called backward here, can still be a forward traversal! See whereToGo() documentation. */
            else if (de_pre->is_not_visited()){
                cout << "\t[ 5] I am at " << de->getID() << " and will go backward to " << de_pre->getID() << endl;
                de->set_visited();
                ret = DFS_Direction_Recursion(predecessor, um, dist);
                cout << "\tI jumped back to ID " << de->getID() << endl;
                continue;
            }
            else{
                cout << "\t[ X] undefined case detected here at ID " << de->getID() << endl;
                return 1;
            }
        }
    }
    else{
        for (auto &successor : um.getSuccessors()){
            DataExtension* de_suc = successor.getData();

            // NOTE: case order matters atm!

            if (successor==src){
                cout << "\t[ 7] small-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (successor==um){
                cout << "\t[ 8] self-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (de_suc->is_visited()){
                cout << "\t[ 9] small bubble detected ending at ID " << de_suc->getID() << endl;
                de->set_visited();
                return Traceback_Init(um, successor);
            }
            else if (de_suc->is_not_visited()){
                cout << "\t[10] I am at " << de->getID() << " and will go forward to " << de_suc->getID() << endl;
                de->set_visited();
                ret = DFS_Direction_Recursion(successor, um, dist);
                cout << "\tI jumped back to ID " << de->getID() << endl;
                continue;
            }
            else{
                cout << "\t[ X] undefined case detected here at ID " << de->getID() << endl;
                return 1;
            }
        }
    }

    // if max dist. was exceeded and mothing more shall happen (valid case)
    return ret;
}


/*!
 * \fn      bool ExtendedCDBG::small_bubble_removal()
 * \brief   This function removes small bubbles (up to length 2k-1) from the graph.
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::small_bubble_removal(){
    cout << "[PROGRESS] Running CDBG small bubble removal..." << endl;

    if (!id_init_status){
        cerr << "ERROR: IDs are not initiated." << endl;
        return 1;
    }

    bool ret = 0;

    for (auto &um : *this){
        ret |= DFS_Direction_Init(um, GO_BACKWARD);     // |= keeps evaluating rhs, || don't
        clear_traversal_attributes();
    }
    for (auto &um : *this){
        ret |= DFS_Direction_Init(um, GO_FORWARD);
        clear_traversal_attributes();
    }

    return ret;
}


/*!
 * \fn      void ExtendedCDBG::clear_traversal_attributes()
 * \brief   This functions resets all traversal variables to their default.
 */
inline void ExtendedCDBG::clear_traversal_attributes(){
    for (auto &um : *this){
        DataExtension* de = um.getData();
        de->set_not_seen_visited();
    }
}


/*!
 * \fn      bool ExtendedCDBG::Traceback_Init(const UnitigMap<DataExtension> &current_um, const UnitigMap<DataExtension> &traceback_src)
 * \brief   This functions initiates a DFS-like traversal to traceback paths in a small bubble found by the directed DFS.
 * \param   traceback_src is the sink of the bubble search
 * \param   src is the unitig from where the bubble property got discovered (direct predecessor of traceback_src).
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::Traceback_Init(const UnitigMap<DataExtension> &src, const UnitigMap<DataExtension> &traceback_src){
    bool ret = 0;

    BubblePathSet bubblePaths;
    BubblePath path;

    uint8_t nb_branchingBubblePaths = 0;

    if (whereToGo(traceback_src, src)==GO_BACKWARD) {       // very important to walk against the logic of whereToGo() for once here: if GO_BACKWARD, traverse successors

        for (auto &successor : traceback_src.getSuccessors()){

            DataExtension* de_suc = successor.getData();
            if (de_suc->is_visited()){
                ret |= Traceback_Visit(successor, traceback_src, path, nb_branchingBubblePaths);
            }

            bubblePaths.push_back(path);
            path.clear();
        }
    }
    else{

        for (auto &predecessor : traceback_src.getPredecessors()){

            DataExtension* de_pre = predecessor.getData();
            if (de_pre->is_visited()){
                ret |= Traceback_Visit(predecessor, traceback_src, path, nb_branchingBubblePaths);
            }

            bubblePaths.push_back(path);
            path.clear();
        }
    }

    if (ret){
        cout << "\tWARNING: A loop occurred during traceback. Algorithm might got confused with directionality." << endl;
        ret = 0;
    }
    else{
        cout << "\t"; prettyprint::print(bubblePaths);
    }

    // TODO: delete a path according to rules

    return ret;
}


bool ExtendedCDBG::Traceback_Visit(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, BubblePath &path, uint8_t &nb_branchingBubblePaths){
    bool ret = 0;

    DataExtension* de = um.getData();
    path.push_back(de->getID());

    if(getDegree(um)>2)
        nb_branchingBubblePaths++;

    if (whereToGo(um, src)==GO_BACKWARD) {
        for (auto &predecessor : um.getPredecessors()){

            if (predecessor==src || predecessor==um)        /* TODO: traceback can't resolve small- or self-loops */
                return 1;

            DataExtension* de_pre = predecessor.getData();
            if (de_pre->is_visited())
                ret |= Traceback_Visit(predecessor, um, path, nb_branchingBubblePaths);
        }

    }
    else{
        for (auto &successor : um.getSuccessors()){

            if (successor==src || successor==um)            /* TODO: traceback can't resolve small- or self-loops */
                return 1;

            DataExtension* de_suc = successor.getData();
            if (de_suc->is_visited())
                ret |= Traceback_Visit(successor, um, path, nb_branchingBubblePaths);
        }

    }

    return ret;
}























