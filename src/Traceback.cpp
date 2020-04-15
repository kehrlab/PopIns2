#include "Traceback.h"



void Traceback::add(const string &unitig, const bool startnode){

    if(!startnode){

        if(_d == 0x1){          // VISIT_PREDECESSOR

            _contig = _contig + unitig.substr(0, unitig.length()-(_k-1));

        }
        else{                   // VISIT_SUCCESSOR (0x0)

            _contig = unitig.substr(_k-1) + _contig;

        }
    }
    else{       // if startnode

        if(_d == 0x1){          // VISIT_PREDECESSOR

            _contig = _contig + unitig.substr(unitig.length()-(_k-1));

        }
        else{                   // VISIT_SUCCESSOR (0x0)

            _contig = unitig.substr(0,_k-1) + _contig;

        }
    }
}


void Traceback::addFullSink(const string &unitig){

    if(_d == 0x1){          // VISIT_PREDECESSOR

        _contig = _contig + unitig;

    }
    else{                   // VISIT_SUCCESSOR (0x0)

        _contig = unitig + _contig;

    }
}


void Traceback::addN(){

    string n_gap(_k, 'N');

    if(_d == 0x1){          // VISIT_PREDECESSOR

        _contig = _contig + n_gap;

    }
    else{                   // VISIT_SUCCESSOR (0x0)

        _contig = n_gap + _contig;

    }
}


void Traceback::write(ofstream &of, const size_t counter) const{

    try{

        of << ">contig_" << counter << "\n";
        of << _contig << "\n";

    }
    catch(std::ofstream::failure &writeErr){

        std::cerr << "\n[Traceback::write] Error:\n"
             << writeErr.what()
             << std::endl;

    }
}
