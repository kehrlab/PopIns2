#include "CompactedDBG_Graph_extension.h"



void ExtendedCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &um : *this){
        DataExtension* de = um.getData();      // de is a POINTER to a UnitigExtension
        de->setID(i);
        ++i;
    }
    init_status = true;
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
     for (auto &predecessor : um.getPredecessors())
         return predecessor == src ? GO_FORWARD : GO_BACKWARD;
}


/*!
 * \fn      bool ExtendedCDBG::BFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction)
 * \brief   This function initiates the recursion of the directed BFS.
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::BFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction){

    bool ret = 0;
    unsigned int dist = 0;

    // mark current unitig
    DataExtension* de = um.getData();
    um.strand==1? cout << "[[" << de->getID() << "+" << "]]" << endl : cout << "[[" << de->getID() << "-" << "]]" << endl;

    if (de->is_not_visited())
        de->set_visited();

    if (direction==GO_BACKWARD){
        for (auto &predecessor : um.getPredecessors()){
            DataExtension* de_pre = predecessor.getData();
            if (predecessor.size < BFS_MAX_DIST){
                cout << "\t[  ] I am at " << de->getID() << " and will go backward to " << de_pre->getID() << endl;
                /*  According to Stack Overflow [https://stackoverflow.com/questions/9021049/operator-for-a-boolean-in-c],
                *  the rhs of |= is always evaluated. That's wanted here since I don't want to break the recursion
                *  in new undefined cases.
                */
                ret |= BFS_Direction_Recursion(predecessor, um, dist);
            }
        }
    }
    else if (direction==GO_FORWARD){
        for (auto &successor : um.getSuccessors()){
            DataExtension* de_suc = successor.getData();
            if (successor.size < BFS_MAX_DIST){
                cout << "\t[  ] I am at " << de->getID() << " and will go forward to " << de_suc->getID() << endl;
                /* see predecessor for equivalent explanation */
                ret |= BFS_Direction_Recursion(successor, um, dist);
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
 * \fn      bool ExtendedCDBG::BFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist)
 * \brief   This function defines the recursion of the directed BFS.
 * \return  bool; 0 for success
 */
bool ExtendedCDBG::BFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist){

    bool ret = 0;
    dist += um.size;

    // unitig is of interest, mark as traversed
    DataExtension* de = um.getData();
    um.strand==1? cout << "\t[[" << de->getID() << "+" << "]]" << endl : cout << "\t[[" << de->getID() << "-" << "]]" << endl;

    if (whereToGo(um, src)==GO_BACKWARD){
        for (auto &predecessor : um.getPredecessors()){
            DataExtension* de_pre = predecessor.getData();

            // NOTE: case order matters atm!

            if (dist>=BFS_MAX_DIST && de_pre->is_not_visited()){
                cout << "\t[ 1] case exceeds max distance here at ID " << de->getID() << endl;
                break;
            }
            else if (predecessor==src){
                cout << "\t[ 2] small-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (predecessor==um){
                cout << "\t[ 3] self-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            /* NOTE: What is called backward here, can still be a forward traversal! See whereToGo() documentation. */
            else if (dist<BFS_MAX_DIST && de_pre->is_not_visited()){
                cout << "\t[ 4] I am at " << de->getID() << " and will go backward to " << de_pre->getID() << endl;
                de->set_visited();
                ret = BFS_Direction_Recursion(predecessor, um, dist);
                cout << "\tI jumped back back to ID " << de->getID() << endl;
                continue;
            }
            else if (de_pre->is_visited()){
                cout << "\t[ 5] small bubble detected ending at ID " << de_pre->getID() << endl;
                de->set_visited();
                // NOTE: DON'T mark sink! In rare cases multiple bubbles might end up in same sink.
                // TODO: trigger traceack
                return 0;
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

            if (dist>=BFS_MAX_DIST && de_suc->is_not_visited()){
                cout << "\t[ 6] case exceeds max distance here at ID " << de->getID() << endl;
                break;
            }
            else if (successor==src){
                cout << "\t[ 7] small-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (successor==um){
                cout << "\t[ 8] self-loop detected here at ID " << de->getID() << endl;
                continue;
            }
            else if (dist<BFS_MAX_DIST && de_suc->is_not_visited()){
                cout << "\t[ 9] I am at " << de->getID() << " and will go forward to " << de_suc->getID() << endl;
                de->set_visited();
                ret = BFS_Direction_Recursion(successor, um, dist);
                cout << "\tI jumped back back to ID " << de->getID() << endl;
                continue;
            }
            else if (de_suc->is_visited()){
                cout << "\t[10] small bubble detected ending at ID " << de_suc->getID() << endl;
                de->set_visited();
                // NOTE: DON'T mark sink! In rare cases multiple bubbles might end up in same sink.
                // TODO: trigger traceack
                return 0;
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
        ret |= BFS_Direction_Init(um, GO_BACKWARD);     // |= keeps evaluating rhs, || don't
        clear_traversal_attributes();
    }
    for (auto &um : *this){
        ret |= BFS_Direction_Init(um, GO_FORWARD);
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



























