#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1


#include <unordered_map>

#include "../src/ColoredDeBruijnGraph.h"
#include "../src/util.h"


CCDBG_Build_opt ccdbg_opt;
SEQAN_DEFINE_TEST(test_ccdbg_opt){

    std::string path = "./simulated/S0001/";
    std::vector<std::string> infiles;
    getFilesFromDir(infiles, path);

    ccdbg_opt.filename_seq_in = infiles;
    ccdbg_opt.deleteIsolated = true;
    ccdbg_opt.clipTips = true;
    ccdbg_opt.prefixFilenameOut = "unit_test";
    ccdbg_opt.nb_threads = 4;
    ccdbg_opt.outputGFA = true;
    ccdbg_opt.verbose = false;
}

ExtendedCCDBG ccdbg(ccdbg_opt.k, ccdbg_opt.g);
SEQAN_DEFINE_TEST(test_ccdbg_build){

    SEQAN_ASSERT_EQ(
        ccdbg.buildGraph(ccdbg_opt),
        true
    );

    SEQAN_ASSERT_EQ(
        ccdbg.simplify(ccdbg_opt.deleteIsolated, ccdbg_opt.clipTips, ccdbg_opt.verbose),
        true
    );

    SEQAN_ASSERT_EQ(
        ccdbg.buildColors(ccdbg_opt),
        true
    );

    /*
    SEQAN_ASSERT_EQ(
        ccdbg.write(ccdbg_opt.prefixFilenameOut, ccdbg_opt.nb_threads, ccdbg_opt.verbose),
        true
    );
    */
}

SEQAN_DEFINE_TEST(test_ccdbg_functions){

    ccdbg.init_ids();
    SEQAN_ASSERT_EQ(
        ccdbg.is_id_init(),
        true
    );

    SEQAN_ASSERT_EQ(
        ccdbg.connected_components(ccdbg_opt),
        true
    );

    SEQAN_ASSERT_EQ(
        ccdbg.count_connected_components(),
        139u
    );
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

    SEQAN_CALL_TEST(test_ccdbg_opt);
    SEQAN_CALL_TEST(test_ccdbg_build);
    SEQAN_CALL_TEST(test_ccdbg_functions);

}
SEQAN_END_TESTSUITE
