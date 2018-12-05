#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1



#include "../src/ColoredDeBruijnGraph.h"
#include "../src/util.h"

using namespace seqan;


CCDBG_Build_opt ccdbg_opt;
SEQAN_DEFINE_TEST(test_ccdbg_opt){

    std::string path = "./simulated/S0001/";
    std::vector<std::string> infiles;
    SEQAN_ASSERT_EQ(getFastx(infiles, path), true);

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

SEQAN_DEFINE_TEST(test_ccdbg_connected_components){

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

CCDBG_Build_opt opt;
ExtendedCCDBG g(opt.k, opt.g);
SEQAN_DEFINE_TEST(test_ccdbg_simpleBranching_singleThread){
    /* Reset graph and input files */
    //ccdbg.clear();
    //ccdbg_opt.filename_seq_in.clear();

    std::string path = "./testcases/simpleBranching/";
    std::vector<std::string> infiles;
    SEQAN_ASSERT_EQ(getFastx(infiles, path), true);

    opt.filename_seq_in = infiles;
    opt.deleteIsolated = true;
    opt.clipTips = true;
    opt.prefixFilenameOut = "simpleBranching";
    opt.nb_threads = 1;
    opt.outputGFA = true;
    opt.verbose = false;

    /* Build and prune graph */
    SEQAN_ASSERT_EQ(g.buildGraph(opt), true);
    SEQAN_ASSERT_EQ(g.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose), true);
    SEQAN_ASSERT_EQ(g.buildColors(opt), true);
    SEQAN_ASSERT_EQ(g.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose), true);

    /* Run merge */
    g.init_ids();
    SEQAN_ASSERT_EQ(g.merge(opt), true);

    /* Truth set */
    StringSet<DnaString> simpleBranchingTruthSet;
    DnaString str1 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCACGATCGCTCTAGCATGCACGGATGTCAGCATGCACATGCGCTTCTTCACGCCCCCCCCA";
    DnaString str2 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGCCATTTAAACCACTCTGGGCCTCAGTTTGCATTAGCCCCCCCTAGTTCGAGCCACACGTGTGTACGTACCGCTAATGCTGGG";
    DnaString str3 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGCCATTTAAACCACTCTGGGCCTTACGCGCGAACTGTACGGGCATAATCGGATCTTTTTCCGATAGTTACCAAACCATGTCGT";
    DnaString str4 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCGTGGACCGAGAGAATACACCACCTGGACCATTGGGCAGTTATTTGAACCAGTCTGACCCTCTACAGTGCTATATATATAACGTAGCGTACGATCATATCGCATCCGTCGCTACGCTATTA";
    appendValue(simpleBranchingTruthSet, str1);
    appendValue(simpleBranchingTruthSet, str2);
    appendValue(simpleBranchingTruthSet, str3);
    appendValue(simpleBranchingTruthSet, str4);

    /* Test contig.fa for correctness */
    CharString seqFileName = "contigs.fa";
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName))) std::cerr << "ERROR: Could not open the file.\n";

    try{
        readRecords(ids, seqs, seqFileIn);
    }
    catch (Exception const & e){
        std::cout << "ERROR: " << e.what() << std::endl;
    }
    SEQAN_ASSERT_EQ(length(seqs), 4u);

    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
        unsigned c = 0;
        for (TStringSetIterator true_seq = begin(simpleBranchingTruthSet); true_seq != end(simpleBranchingTruthSet); ++true_seq)
            if (*seq == *true_seq)
                ++c;
        SEQAN_ASSERT_EQ(c, 1u);
    }
    
    if(remove("contigs.fa") || remove("simpleBranching.gfa") || remove("simpleBranching.bfg_colors"))
        perror("Error deleting file");
    else
        puts("Files successfully deleted");
}

CCDBG_Build_opt opt2;
ExtendedCCDBG g2(opt2.k, opt2.g);
SEQAN_DEFINE_TEST(test_ccdbg_simpleBubbles_singleThread){
    /* Reset graph and input files */
    //ccdbg.clear();
    //ccdbg_opt.filename_seq_in.clear();

    std::string path = "./testcases/simpleBubbles/";
    std::vector<std::string> infiles;
    SEQAN_ASSERT_EQ(getFastx(infiles, path), true);

    opt2.filename_seq_in = infiles;
    opt2.deleteIsolated = true;
    opt2.clipTips = true;
    opt2.prefixFilenameOut = "simpleBubbles";
    opt2.nb_threads = 1;
    opt2.outputGFA = true;
    opt2.verbose = false;

    /* Build and prune graph */
    SEQAN_ASSERT_EQ(g2.buildGraph(opt2), true);
    SEQAN_ASSERT_EQ(g2.simplify(opt2.deleteIsolated, opt2.clipTips, opt2.verbose), true);
    SEQAN_ASSERT_EQ(g2.buildColors(opt2), true);
    SEQAN_ASSERT_EQ(g2.write(opt2.prefixFilenameOut, opt2.nb_threads, opt2.verbose), true);

    /* Run merge */
    g2.init_ids();
    SEQAN_ASSERT_EQ(g2.merge(opt2), true);

    /* Truth set */
    StringSet<DnaString> simpleBubblesTruthSet;
    DnaString str1 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGTTATTTAAACCAGTCTGAGCCTCCCCCAGAGCCGTTCGCGCCGCCCCCGGTCCTCCGGCCCCCGGTCTGCCCCGCAGCGCCTGCCCGGCTTAATGTCAGAGACAGCCCACCCACTCCATAAATCCACTTGTGACAGGGCTGGGGACCTGGACTGTCCTCAGAGAGGCCCCCTGTGACCACTC";
    DnaString str2 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGTTATTTAAACCAGTCTGAGCCTCCCCCAGAGCCGTTCGCGCCGCCCCCGGTCCTCCGGCCCCCGGTCTGCCCCGCAGCGCCTGCCCGGCTTGTTTTGGTATTATTCATCGTGAGGTGAAGACCAAATTTCTCCTCAGAGATGCAAGGGCTACGT";
    appendValue(simpleBubblesTruthSet, str1);
    appendValue(simpleBubblesTruthSet, str2);

    /* Test contig.fa for correctness */
    CharString seqFileName = "contigs.fa";
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName))) std::cerr << "ERROR: Could not open the file.\n";

    try{
        readRecords(ids, seqs, seqFileIn);
    }
    catch (Exception const & e){
        std::cout << "ERROR: " << e.what() << std::endl;
    }
    SEQAN_ASSERT_EQ(length(seqs), 2u);

    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
        unsigned c = 0;
        for (TStringSetIterator true_seq = begin(simpleBubblesTruthSet); true_seq != end(simpleBubblesTruthSet); ++true_seq)
            if (*seq == *true_seq)
                ++c;
        SEQAN_ASSERT_EQ(c, 1u);
    }
    
    if(remove("contigs.fa") || remove("simpleBubbles.gfa") || remove("simpleBubbles.bfg_colors"))
        perror("Error deleting file");
    else
        puts("Files successfully deleted");
}

SEQAN_BEGIN_TESTSUITE(test_popins2){

    SEQAN_CALL_TEST(test_ccdbg_opt);
    SEQAN_CALL_TEST(test_ccdbg_build);
    SEQAN_CALL_TEST(test_ccdbg_connected_components);
    SEQAN_CALL_TEST(test_ccdbg_simpleBranching_singleThread);
    SEQAN_CALL_TEST(test_ccdbg_simpleBubbles_singleThread);

}
SEQAN_END_TESTSUITE
