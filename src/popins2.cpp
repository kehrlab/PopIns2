#include <iostream>
#include <bifrost/CompactedDBG.hpp>
#include "argument_parsing.h"           /* seqAn argument parser */
#include "CDBG_Data_extension.h"
#include "ExtendedCDBG.h"
using namespace std;


static const uint8_t bitmask_encoder_successor[256] = {
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00001100, 0b00000000, 0b00001101, 0b00000000, 0b00000000, 0b00000000, 0b00001110, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,     // t[4][1]='A', t[4][3]='C', t[4][7]='G'
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00001111, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,     // t[5][4]='T'
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000
};


static const uint8_t bitmask_encoder_predecessor[256] = {
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000011, 0b00000000, 0b00000111, 0b00000000, 0b00000000, 0b00000000, 0b00001011, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,     // t[4][1]='A', t[4][3]='C', t[4][7]='G'
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00001111, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,     // t[5][4]='T'
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
  0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000
};


static const char bitmask_decoder[16] = {
    0b01001110, 0b01001110, 0b01001110, 0b01000001,     // N, N, N, A
    0b01001110, 0b01001110, 0b01001110, 0b01000011,     // N, N, N, C
    0b01001110, 0b01001110, 0b01001110, 0b01000111,     // N, N, N, G
    0b01000001, 0b01000011, 0b01000111, 0b01010100      // A, C, G, T
};


// ==============================
// Function: main()
// ==============================
int main(int argc, char const *argv[]){

    // ==============================
    // Argument Parser
    // ==============================
    OptionsWrapper options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // catch parse error
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // ==============================
    // Global program variables
    // ==============================
    // all file names in --indir with full path
    vector<string> sample_fastx_names;

    // ==============================
    // Loop for all FastX files
    // ==============================
    bool isFilereadSuccessful = detect_indir_files(options, sample_fastx_names);

    // ==============================
    // Initialize graph
    // ==============================
    CDBG_Build_opt graph_options;
    if(isFilereadSuccessful==EXIT_SUCCESS)
        init_graph_options(options, sample_fastx_names, graph_options);
    if(!check_ProgramOptions(graph_options))
        return 1;   // some input options are not appropriate to construct the CDBG

    // ==============================
    // Run graph functions
    // ==============================
    ExtendedCDBG cdbg(graph_options.k, graph_options.g);
    cout << "[PROGRESS] Building CDBG..." << endl;
    cdbg.build(graph_options);

    // simplify
    cout << "[PROGRESS] Simplifying CDBG..." << endl;
    cdbg.simplify(graph_options.deleteIsolated, graph_options.clipTips, graph_options.verbose);

    // TEST START
    cdbg.init_ids();

    uint8_t test_bit_encoding;
    for (auto &unitig : cdbg){
        ForwardCDBG<UnitigExtension, false> fw_dbg = unitig.getSuccessors();
        BackwardCDBG<UnitigExtension, false> bw_dbg = unitig.getPredecessors();

        unsigned i=0;
        for (auto &suc : fw_dbg) ++i;
        unsigned j=0;
        for (auto &pre : bw_dbg) ++j;

        if (unitig.getData()->getID() == 4) {
            UnitigMap<UnitigExtension>::neighbor_iterator suc1 = fw_dbg.begin();
            UnitigMap<UnitigExtension>::neighbor_iterator pre1 = bw_dbg.begin();

            size_t pre1id = pre1->getData()->getID();
            cout << "Predecessor: " << pre1id << endl;
            char pre1tail = (pre1->getTail()).getChar(graph_options.k-1);   // getChar(offset) returns the base from Kmer at offset from the beginning
            cout << "Predecessor's last base: " << pre1tail << endl;

            size_t suc1id = suc1->getData()->getID();
            cout << "Successor: " << suc1id << endl;
            char suc1head = (suc1->getHead()).getChar(0);                   // getChar(offset) returns the base from Kmer at offset from the beginning
            cout << "Successor's first base: " << suc1head << endl;

            uint8_t pre_bit = bitmask_encoder_predecessor[pre1tail];
            uint8_t suc_bit = bitmask_encoder_successor[suc1head];
            test_bit_encoding = setNeighborPairFromBases(pre_bit, suc_bit);
        }
    }
    pair<uint8_t, uint8_t> test_bit_decoding = getBasesFromNeighborPair(test_bit_encoding);
    uint8_t f = test_bit_decoding.first;
    uint8_t s = test_bit_decoding.second;
    cout << "\n After bit operations the PRE and SUC are: (" << bitmask_decoder[f] << "," << bitmask_decoder[s] << ")" << endl;

    // TEST END

    // write
    cout << "[PROGRESS] Writing GFA..." << endl;
    cdbg.write(graph_options.prefixFilenameOut, graph_options.nb_threads, true, graph_options.verbose);
    cout << "[DEBUG] The DBG has " << cdbg.size() << " unitigs.\n" << endl;

    return 0;
}



