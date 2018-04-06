#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1


#include <unordered_map>

#include "../src/ColoredCDBG_Graph_extension.h"
#include "../src/argument_parsing.h"



CCDBG_Build_opt graph_opt;
SEQAN_DEFINE_TEST(init_bifrost_parameter){
    // ---------------------------------
    // | Initiate graph build options  |
    // ---------------------------------
    std::vector<std::string> infiles;
    infiles.push_back("./unitigs/simulated_S0001.fasta");
    graph_opt.filename_seq_in = infiles;
    graph_opt.prefixFilenameOut = "union_test_out";
    graph_opt.nb_threads = 4;
}

ExtendedCCDBG xg(graph_opt.k, graph_opt.g);
SEQAN_DEFINE_TEST(test_bifrost_graphfunctions){
    // -----------------------
    // | Test Bifrost methods  |
    // -----------------------
    xg.build(graph_opt);
    SEQAN_ASSERT_EQ(xg.size(), 309u);
    xg.mapColors(graph_opt);
    SEQAN_ASSERT_EQ(xg.write(graph_opt.prefixFilenameOut, graph_opt.nb_threads, graph_opt.verbose), true);
}

SEQAN_DEFINE_TEST(test_init_ids){
    xg.init_ids();
    SEQAN_ASSERT_EQ(xg.is_init(), true);
}

SEQAN_DEFINE_TEST(test_connectedcomponents){
    bool cc_build = xg.connected_components(graph_opt);
    SEQAN_ASSERT_EQ(cc_build, true);
    size_t nb_cc = xg.count_connected_components();
    SEQAN_ASSERT_EQ(nb_cc, 139u);
}

// NOTE: TEST was designed for other testfile. DO NOT RUN!
/*
SEQAN_DEFINE_TEST(test_neighbors_and_bit_operations){
    size_t nb_potential_splitnodes = 0;
    for (auto &ucm : xg){
        // -----------------------------
        // | Test amount of neighbors  |
        // -----------------------------
        ForwardCDBG<UnitigExtension, false> fw_dbg = ucm.getSuccessors();
        BackwardCDBG<UnitigExtension, false> bw_dbg = ucm.getPredecessors();

        size_t i=0;
        for (auto &suc : fw_dbg) ++i;

        size_t j=0;
        for (auto &pre : bw_dbg) ++j;

        if (i>1 && j>1)
            ++nb_potential_splitnodes;

        // ------------------------
        // | Test base accession  |
        // ------------------------
        if (ucm.getData()->getID() == 4) {
            UnitigMap<UnitigExtension>::neighbor_iterator suc1 = fw_dbg.begin();
            UnitigMap<UnitigExtension>::neighbor_iterator pre1 = bw_dbg.begin();

            char firstPredecessor_lastBase = (pre1->getTail()).getChar(graph_opt.k-1);   // getChar(offset) returns the base from Kmer at offset from the beginning
            char firstSuccessor_firstBase = (suc1->getHead()).getChar(0);                   // getChar(offset) returns the base from Kmer at offset from the beginning

            SEQAN_ASSERT_EQ(firstPredecessor_lastBase, 'T');
            SEQAN_ASSERT_EQ(firstSuccessor_firstBase, 'C');

            // -------------------------------------------
            // | Test succint bit storage and bit masks  |
            // -------------------------------------------
            uint8_t pre_bitmask = bitmask_encoder_predecessor[firstPredecessor_lastBase];
            uint8_t suc_bitmask = bitmask_encoder_successor[firstSuccessor_firstBase];
            uint8_t test_bit_encoding = setNeighborPairFromBases(pre_bitmask, suc_bitmask);
            SEQAN_ASSERT_EQ(test_bit_encoding, 0b00001101);

            pair<uint8_t, uint8_t> test_bit_decoding = getBasesFromNeighborPair(test_bit_encoding);
            char f = bitmask_decoder[test_bit_decoding.first];
            char s = bitmask_decoder[test_bit_decoding.second];
            SEQAN_ASSERT_EQ(f, 'T');
            SEQAN_ASSERT_EQ(s, 'C');
        }
    }
    SEQAN_ASSERT_EQ(nb_potential_splitnodes, 5u);
}
*/


SEQAN_BEGIN_TESTSUITE(test_popins2){

	// call tests here
    SEQAN_CALL_TEST(init_bifrost_parameter);
    SEQAN_CALL_TEST(test_bifrost_graphfunctions);
    SEQAN_CALL_TEST(test_init_ids);
    SEQAN_CALL_TEST(test_connectedcomponents);
    //SEQAN_CALL_TEST(test_neighbors_and_bit_operations);

}
SEQAN_END_TESTSUITE
