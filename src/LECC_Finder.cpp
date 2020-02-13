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

    return LECC_;
}


void LECC_Finder::annotate_recursion(UnitigColorMap<UnitigExtension> &ucm, const unsigned LECC__){
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    // skip if already LECC-annotated (0 is default LECC identifier)
    if(ue->getLECC() != 0)
        return;

    // skip high entropy unitigs
    if(ue->getEntropy() >= this->threshold_)
        return;

    // found new LECC
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
