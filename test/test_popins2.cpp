#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1



#include "../src/ColoredDeBruijnGraph.h"
#include "../src/util.h"

using namespace seqan;

typedef std::vector<std::string> strings_v;
typedef uint8_t direction_t;

const static direction_t VISIT_SUCCESSOR   = 0x0;
const static direction_t VISIT_PREDECESSOR = 0x1;

template <typename TType>
inline void print(std::vector<TType> &v){
    std::cout << "[";
    // stay safe with constant iterators
    typename std::vector<TType>::const_iterator it = v.cbegin();
    for ( ; it != v.cend(); ++it){
        // output format
        std::cout << *it << ", ";
    }
    std::cout << "]" << std::endl;
}


/** Get the color overlap of two neighbor unitig.
* The function isolates the kmers that face each other with respect to the unititgs. Then, it retrieves
* the color vectors of both kmers,does an AND operation and counts the intersecton.
* @param ucm_to_get_head_from is the unitig to get the head kmer from
* @param ucm_to_get_tail_from is the unitig to get the tail kmer from
* @param nb_colors is the number of samples in the graph (NOTE: might be obsolete in non-test code, use CCDBG->getNbColors())
**/
inline unsigned get_neighbor_overlap(const UnitigColorMap<UnitigExtension> &ucm_to_get_head_from, const UnitigColorMap<UnitigExtension> &ucm_to_get_tail_from, const unsigned nb_colors){
        size_t len = ucm_to_get_tail_from.len;   // I assume this gets me the past-last-kmer index
        //std::cout << len << std::endl;

        const UnitigColorMap<UnitigExtension> k_first = ucm_to_get_head_from.getKmerMapping(0);
        const UnitigColorMap<UnitigExtension> k_last  = ucm_to_get_tail_from.getKmerMapping(len-1);

        const UnitigColors* k_first_colors = k_first.getData()->getUnitigColors(k_first);
        const UnitigColors* k_last_colors  =  k_last.getData()->getUnitigColors(k_last);

        std::vector<bool> k_first_color_bits(nb_colors, false);
        std::vector<bool> k_last_color_bits(nb_colors, false);

        // get color IDs of unitig's first kmer
        UnitigColors::const_iterator cit = k_first_colors->begin(k_first);
        for (; cit != k_first_colors->end(); ++cit)
            k_first_color_bits[cit.getColorID()] = true;

        // get color IDs of unitig's last kmer
        cit = k_last_colors->begin(k_last);
        for (; cit != k_last_colors->end(); ++cit)
            k_last_color_bits[cit.getColorID()] = true;

        // PRINT
        //std::string first_seq = k_first.mappedSequenceToString();
        //std::string last_seq  =  k_last.mappedSequenceToString();
        //std::cout << "START COLORS (" << first_seq << "): " ; print(k_first_color_ids);
        //std::cout << "END COLORS ("   << last_seq  << "): " ; print(k_last_color_ids);

        // sum of intersection
        unsigned count = 0;
        for (unsigned i=0; i < k_first_color_bits.size(); ++i)                  // NOTE: IMPROVEMENT: comparisons can be reduced by using UnitigColors::colorMax(ucm)
            if (k_first_color_bits[i] && k_last_color_bits[i])
                ++count;

        return count;
}


/** Get the ID of the best fitting neighbor.
* This function iterates over all neighbors with respect to the traversal direction. It then applies
* the function get_neighbor_overlap() to determine the neighbor that has the highest color match.
* @param ucm ist the unitig to compare to
* @param neighbors is a ForwardCDBG or BackwardCDBG, its unitigs will be compared to ucm
* @param direction is the traversal direction
* @param nb_colors is the number of samples in the graph (NOTE: might be obsolete in non-test code, use CCDBG->getNbColors())
**/
template <typename TNeighbors>
inline unsigned get_best_neighbor_ID(const UnitigColorMap<UnitigExtension> &ucm, const TNeighbors &neighbors, const direction_t direction, const unsigned nb_colors){

    unsigned best_neighbor_id = 0; // default
    unsigned best_overlap = 0;

    if (direction==VISIT_PREDECESSOR){

        for (auto &pre : neighbors){

            unsigned overlap = get_neighbor_overlap(ucm, pre, nb_colors);

            if (overlap > best_overlap){

                DataAccessor<UnitigExtension>* da = pre.getData();              // NOTE: can I avoid recreating these pointer over and over again?
                UnitigExtension* data = da->getData(pre);
                best_neighbor_id = data->getID();
            }
        }
    }
    else{   // direction==VISIT_SUCCESSOR

        for (auto &suc : neighbors){

            unsigned overlap = get_neighbor_overlap(suc, ucm, nb_colors);

            if (overlap > best_overlap){

                DataAccessor<UnitigExtension>* da = suc.getData();              // NOTE: can I avoid recreating these pointer over and over again?
                UnitigExtension* data = da->getData(suc);
                best_neighbor_id = data->getID();
            }
        }
    }

    return best_neighbor_id;
}


inline bool is_startnode(const UnitigColorMap<UnitigExtension> &ucm){
    return
        ucm.len>2 &&                                                                                    // be longer than 2 kmers and
        (( ucm.getPredecessors().hasPredecessors() && !ucm.getSuccessors().hasSuccessors() ) ||         // have only predecessors or
         (!ucm.getPredecessors().hasPredecessors() &&  ucm.getSuccessors().hasSuccessors() ) );         // have only successors
}


/**
* @return 1 for further traversal, 0 for jump back into parent recursion level
**/
inline uint8_t DFS(UnitigColorMap<UnitigExtension> &ucm, const direction_t direction){

    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* data = da->getData(ucm);

    if(direction==VISIT_PREDECESSOR){

        if(data->is_undiscovered_bw()){

            data->set_seen_bw();

            BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> predecessors = ucm.getPredecessors();

            if (!predecessors.hasPredecessors()){       // sink node

                return 0;
            }
        }
    }
    else{   // if direction==VISIT_SUCCESSOR

        if(data->is_undiscovered_fw()){

            data->set_seen_fw();

            ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> successors = ucm.getSuccessors();

            if (!successors.hasSuccessors()){       // sink node

                return 0;
            }
        }
    }
}


/** This function traverses the graph.
* @return 1 for successful execution
**/
inline uint8_t traverse(ExtendedCCDBG &g){

    if (!g.is_id_init()) return 0;      // sanity check

    for (auto &ucm : g){

        if (is_startnode(ucm)){             // NOTE: improvement: receive traversal direction from is_startnode() s.t. traverse() doesn't have to find out again

            if (ucm.getPredecessors().hasPredecessors()){

                DFS(ucm, VISIT_PREDECESSOR);
            }
            else{ //ucm.getSuccessors().hasSuccessors()

                DFS(ucm, VISIT_SUCCESSOR);
            }
        }
    }
}


// -----------------
// | LECC UNITTEST |
// -----------------
CCDBG_Build_opt opt_lecc_unittest;

SEQAN_DEFINE_TEST(setup_lecc_unittest){
    strings_v i_files;
    getFastx(i_files, "./testcases/lecc_unittest/");

    opt_lecc_unittest.filename_ref_in = i_files;
    opt_lecc_unittest.deleteIsolated = true;
    opt_lecc_unittest.clipTips = true;
    opt_lecc_unittest.prefixFilenameOut = "lecc_unittest";
    opt_lecc_unittest.nb_threads = 1;
    opt_lecc_unittest.outputGFA = true;
    opt_lecc_unittest.verbose = false;
    opt_lecc_unittest.k = 63;
}

ExtendedCCDBG ccdbg_lecc_unittest(opt_lecc_unittest.k, opt_lecc_unittest.g);

SEQAN_DEFINE_TEST(test_lecc_unittest){
    SEQAN_ASSERT_EQ(ccdbg_lecc_unittest.buildGraph(opt_lecc_unittest),true);

    SEQAN_ASSERT_EQ(ccdbg_lecc_unittest.simplify(opt_lecc_unittest.deleteIsolated, opt_lecc_unittest.clipTips, opt_lecc_unittest.verbose),true);

    SEQAN_ASSERT_EQ(ccdbg_lecc_unittest.buildColors(opt_lecc_unittest),true);



    SEQAN_ASSERT_EQ(ccdbg_lecc_unittest.write(opt_lecc_unittest.prefixFilenameOut, opt_lecc_unittest.nb_threads, opt_lecc_unittest.verbose), true);
}







// =========== legacy content ==============
/*
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

// -------------
// | NEW GRAPH |
// -------------

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


    //SEQAN_ASSERT_EQ(
    //    ccdbg.write(ccdbg_opt.prefixFilenameOut, ccdbg_opt.nb_threads, ccdbg_opt.verbose),
    //    true
    //);

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

// -------------
// | NEW GRAPH |
// -------------
//CCDBG_Build_opt opt;
//ExtendedCCDBG g(opt.k, opt.g);
SEQAN_DEFINE_TEST(test_ccdbg_simpleBranching_singleThread){
    // Reset graph and input files
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

    // Build and prune graph
    SEQAN_ASSERT_EQ(g.buildGraph(opt), true);
    SEQAN_ASSERT_EQ(g.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose), true);
    SEQAN_ASSERT_EQ(g.buildColors(opt), true);
    SEQAN_ASSERT_EQ(g.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose), true);

    // Run merge
    g.init_ids();
    SEQAN_ASSERT_EQ(g.merge(opt), true);

    // Truth set
    StringSet<DnaString> simpleBranchingTruthSet;
    DnaString str1 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCACGATCGCTCTAGCATGCACGGATGTCAGCATGCACATGCGCTTCTTCACGCCCCCCCCA";
    DnaString str2 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGCCATTTAAACCACTCTGGGCCTCAGTTTGCATTAGCCCCCCCTAGTTCGAGCCACACGTGTGTACGTACCGCTAATGCTGGG";
    DnaString str3 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGCCATTTAAACCACTCTGGGCCTTACGCGCGAACTGTACGGGCATAATCGGATCTTTTTCCGATAGTTACCAAACCATGTCGT";
    DnaString str4 = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCGTGGACCGAGAGAATACACCACCTGGACCATTGGGCAGTTATTTGAACCAGTCTGACCCTCTACAGTGCTATATATATAACGTAGCGTACGATCATATCGCATCCGTCGCTACGCTATTA";
    appendValue(simpleBranchingTruthSet, str1);
    appendValue(simpleBranchingTruthSet, str2);
    appendValue(simpleBranchingTruthSet, str3);
    appendValue(simpleBranchingTruthSet, str4);

    // Test contig.fa for correctness
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
        perror("[test_ccdbg_simpleBranching_singleThread] Error deleting meta files");
    //else
    //    puts("[test_ccdbg_simpleBranching_singleThread] Meta files successfully deleted");
}

// -------------
// | NEW GRAPH |
// -------------

CCDBG_Build_opt opt2;
ExtendedCCDBG g2(opt2.k, opt2.g);
SEQAN_DEFINE_TEST(test_ccdbg_simpleBubbles){
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

    // Build and prune graph
    SEQAN_ASSERT_EQ(g2.buildGraph(opt2), true);
    SEQAN_ASSERT_EQ(g2.simplify(opt2.deleteIsolated, opt2.clipTips, opt2.verbose), true);
    SEQAN_ASSERT_EQ(g2.buildColors(opt2), true);
    SEQAN_ASSERT_EQ(g2.write(opt2.prefixFilenameOut, opt2.nb_threads, opt2.verbose), true);

    // Run merge
    g2.init_ids();
    SEQAN_ASSERT_EQ(g2.merge(opt2), true);

    // Truth set
    StringSet<DnaString> simpleBubblesTruthSet;
    DnaString longPathWithTwoBubbles1   = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCCTGGAAAGAGAGAAGACACCAGCTGGACAATTCAGCAGTTATTTAAACCAGTCTGAGCCTCCCCCAGAGCCGTTCGCGCCGCCCCCGGTCCTCCGGCCCCCGGTCTGCCCCGCAGCGCCTGCCCGGCTTAATGTCAGAGACAGCCCACCCACTCCATAAATCCACTTGTGACAGGGCTGGGGACCTGGACTGTCCTCAGAGAGGCCCCCTGTGACCACTC";
    DnaString longPathWithTwoBubbles1rc = "GAGTGGTCACAGGGGGCCTCTCTGAGGACAGTCCAGGTCCCCAGCCCTGTCACAAGTGGATTTATGGAGTGGGTGGGCTGTCTCTGACATTAAGCCGGGCAGGCGCTGCGGGGCAGACCGGGGGCCGGAGGACCGGGGGCGGCGCGAACGGCTCTGGGGGAGGCTCAGACTGGTTTAAATAACTGCTGAATTGTCCAGCTGGTGTCTTCTCTCTTTCCAGGATCAGGAGTTTGAGACCAGCCTGGCTAACATGATGAAACCCTGTCTTTACTAAAAATACAATGATTAGCTGGGTGTGGTGGCGGGCGCCTGTAATTCCAG";
    DnaString longPathWithTwoBubbles2   = "CTGGAATTACAGGCGCCCGCCACCACACCCAGCTAATCATTGTATTTTTAGTAAAGACAGGGTTTCATCATGTTAGCCAGGCTGGTCTCAAACTCCTGATCGTGGACCGAGAGAATACACCACCTGGACAATTCAGCAGTTATTTAAACCAGTCTGAGCCTCCCCCAGAGCCGTTCGCGCCGCCCCCGGTCCTCCGGCCCCCGGTCTGCCCCGCAGCGCCTGCCCGGCTTCTCCGACTCAAAGTGACCAGCACCTCATAAATCCACTTGTGACAGGGCTGGGGACCTGGACTGTCCTCAGAGAGGCCCCCTGTGACCACTC";
    DnaString longPathWithTwoBubbles2rc = "GAGTGGTCACAGGGGGCCTCTCTGAGGACAGTCCAGGTCCCCAGCCCTGTCACAAGTGGATTTATGAGGTGCTGGTCACTTTGAGTCGGAGAAGCCGGGCAGGCGCTGCGGGGCAGACCGGGGGCCGGAGGACCGGGGGCGGCGCGAACGGCTCTGGGGGAGGCTCAGACTGGTTTAAATAACTGCTGAATTGTCCAGGTGGTGTATTCTCTCGGTCCACGATCAGGAGTTTGAGACCAGCCTGGCTAACATGATGAAACCCTGTCTTTACTAAAAATACAATGATTAGCTGGGTGTGGTGGCGGGCGCCTGTAATTCCAG";
    appendValue(simpleBubblesTruthSet, longPathWithTwoBubbles1);
    appendValue(simpleBubblesTruthSet, longPathWithTwoBubbles1rc);
    appendValue(simpleBubblesTruthSet, longPathWithTwoBubbles2);
    appendValue(simpleBubblesTruthSet, longPathWithTwoBubbles2rc);

    // Test contig.fa for correctness //
    CharString seqFileName = "contigs.fa";      // intermediate result file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName))) std::cerr << "ERROR: Could not open the file.\n";
    try{readRecords(ids, seqs, seqFileIn);}
    catch (Exception const & e){std::cout << "ERROR: " << e.what() << std::endl;}
    SEQAN_ASSERT_EQ(length(seqs), 2u);      // intermediate results have to contain two assembled strings
    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
        unsigned c = 0;
        for (TStringSetIterator true_seq = begin(simpleBubblesTruthSet); true_seq != end(simpleBubblesTruthSet); ++true_seq)
            if (*seq == *true_seq)
                ++c;
        SEQAN_ASSERT_EQ(c, 1u);
    }

    //Clean up//
    if(remove("contigs.fa") || remove("contigs.csv") || remove("simpleBubbles.gfa") || remove("simpleBubbles.bfg_colors"))
        perror("[test_ccdbg_simpleBubbles] Error deleting meta files");
}

// -------------
// | NEW GRAPH |
// -------------

CCDBG_Build_opt opt3;
ExtendedCCDBG g3(opt3.k, opt3.g);
SEQAN_DEFINE_TEST(test_ccdbg_simpleColorTest){

    std::string path = "./testcases/simpleColorTest/";
    std::vector<std::string> infiles;
    SEQAN_ASSERT_EQ(getFastx(infiles, path), true);

    opt3.filename_seq_in = infiles;
    opt3.deleteIsolated = true;
    opt3.clipTips = true;
    opt3.prefixFilenameOut = "simpleColorTest";
    opt3.nb_threads = 1;
    opt3.outputGFA = true;
    opt3.verbose = false;

    // Build and prune graph
    SEQAN_ASSERT_EQ(g3.buildGraph(opt3), true);
    SEQAN_ASSERT_EQ(g3.simplify(opt3.deleteIsolated, opt3.clipTips, opt3.verbose), true);
    SEQAN_ASSERT_EQ(g3.buildColors(opt3), true);
    SEQAN_ASSERT_EQ(g3.write(opt3.prefixFilenameOut, opt3.nb_threads, opt3.verbose), true);

    //Run merge //
    g3.init_ids();
    SEQAN_ASSERT_EQ(g3.merge(opt3), true);

    // Truth set //
    StringSet<DnaString> simpleColorTestTruthSet;
    DnaString str1 = "TTACCTCTACAAAAAGCGACTGCCAGTGTAACCCCACGAGGATCCGAAAAGGCGAACCGGCCCAGACGACCCGGGGGCACGGGCCTCAAAGCCGCGACACGACGGCTGTCGGCCGGTAACAGTAACCCCGGAGTGAACTCCTATGGGGCTGGATAGAACAGCCCTGGTGGGCCCCATCAGCAACCCGAATACGTGGCTCTTCGGGAGGCGGCCGGAGGGGCGATGTCTTCCACTATTCGAGGCCGTTCGTTAATACTTGTTGCGTTCCTAGCCGCTATATTTGTCTCTTTGCCGACTAATGTGGACAAGCACACCATAGCCATTTGTCGGGGCGCCTCGGAATACGGTATGAGCAGGCGCCTCGTGAGGCCATCGCGAATACCAGGTGTCCTGTAAGCAGCGAAGGCCCGCACGCGAGATAAACTGCTAGGGAACCGCGTGTCCACGACCGGTGGTGGATTTAATCTCGCCGACGTGTAGACATTCCAGGCAGTGCGTCT";
    DnaString str2 = "GCCGCCGGGCCCCTCTGGTGACTGGGTAGCTGGACTTGCCCTTGGAAGACATAGCAAGACCCTGCCTCTCTATTGATGTCACGGCGAATGTCGGGGAGACAGCAGCGGCTGCAGACATCAGATAACCCCGGAGTGAACTCCTATGGGGCTGGATAGAACAGCCCTGGTGGGCCCCATCAGCAACCCGAATACGTGGCTCTTCGGGAGGCGGCCGGAGGGGCGATGTCTTCCACTATTCGAGGCCGTTCGTTAATACTTGTTGCGTTCCTAGCCGCTATATTTGTCTCTTTGCCGACTAATGTGGACAAGCACACCATAGCCATTTGTCGGGGCGCCTCGGAATACGGTATGAGCAGGCGCCTCGTGCGTTACCGCCATAAGATGGGAGCATGACTCCTTCTCCGCTGCGCCCACGCCAGTAGTGATTACTCCTATGACCCTTCCGAGAGTCCGGAGGCGGAAATCCGCCACGAATGAGAATGTATTTCCCCGACAGTCAT";
    appendValue(simpleColorTestTruthSet, str1);
    appendValue(simpleColorTestTruthSet, str2);

    // Test contig.fa for correctness //
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
    SEQAN_ASSERT_EQ(length(ids), 2u);

    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
        unsigned c = 0;
        for (TStringSetIterator true_seq = begin(simpleColorTestTruthSet); true_seq != end(simpleColorTestTruthSet); ++true_seq)
            if (*seq == *true_seq)
                ++c;
        SEQAN_ASSERT_EQ(c, 1u);
    }

    if(remove("contigs.fa") || remove("simpleColorTest.gfa") || remove("simpleColorTest.bfg_colors"))
        perror("[test_ccdbg_simpleColorTest] Error deleting meta files");
    //else
    //    puts("[test_ccdbg_simpleColorTest] Meta files successfully deleted");
}


CCDBG_Build_opt opt4;
ExtendedCCDBG g4(opt4.k, opt4.g);
SEQAN_DEFINE_TEST(test_ccdbg_rev_comp){

    // *********************************************
    // * This test only checks if Bifrost is using *
    // * the right orientation for unitigs, not if *
    // * the traversal is concatenating correctly. *
    // *********************************************

    std::string path = "./testcases/simpleReverseComplementAssembly/";
    std::vector<std::string> infiles;
    SEQAN_ASSERT_EQ(getFastx(infiles, path), true);

    opt4.filename_seq_in = infiles;
    opt4.deleteIsolated = true;
    opt4.clipTips = true;
    opt4.prefixFilenameOut = "simpleRevComp";
    opt4.nb_threads = 1;
    opt4.outputGFA = true;
    opt4.verbose = true;

    // Build and prune graph //
    SEQAN_ASSERT_EQ(g4.buildGraph(opt4), true);
    SEQAN_ASSERT_EQ(g4.simplify(opt4.deleteIsolated, opt4.clipTips, opt4.verbose), true);
    SEQAN_ASSERT_EQ(g4.buildColors(opt4), true);
    SEQAN_ASSERT_EQ(g4.write(opt4.prefixFilenameOut, opt4.nb_threads, opt4.verbose), true);

    // Run merge //
    g4.init_ids();
    SEQAN_ASSERT_EQ(g4.merge(opt4), true);

    // Truth set //
    StringSet<DnaString> simpleRevCompTruthSet;
    DnaString trueAssembly    = "GAGATGGCAGCGACTAACTCACAGCGCACGTTTGGCTGACGCATATGCATTAGAATCGTCGACGAGCGAGATAGTCTCCGCTGTCGAGGGACTTGCAGAC";
    DnaString trueAssembly_rc = "GTCTGCAAGTCCCTCGACAGCGGAGACTATCTCGCTCGTCGACGATTCTAATGCATATGCGTCAGCCAAACGTGCGCTGTGAGTTAGTCGCTGCCATCTC";
    appendValue(simpleRevCompTruthSet, trueAssembly);
    appendValue(simpleRevCompTruthSet, trueAssembly_rc);

    // Test contig.fa for correctness //
    CharString seqFileName = "contigs.fa";      // intermediate result file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName))) std::cerr << "ERROR: Could not open the file.\n";
    try{readRecords(ids, seqs, seqFileIn);}
    catch (Exception const & e){std::cout << "ERROR: " << e.what() << std::endl;}

    SEQAN_ASSERT_EQ(length(seqs), 1u);      // intermediate result has to contain one assembled string

    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
        unsigned c = 0;
        for (TStringSetIterator true_seq = begin(simpleRevCompTruthSet); true_seq != end(simpleRevCompTruthSet); ++true_seq)
            if (*seq == *true_seq)
                ++c;
        SEQAN_ASSERT_EQ(c, 1u);
    }

    //Clean up//
    if(remove("contigs.fa") || remove("contigs.csv") || remove("simpleRevComp.gfa") || remove("simpleRevComp.bfg_colors"))
        perror("[test_ccdbg_rev_comp] Error deleting meta files");
}
*/


SEQAN_BEGIN_TESTSUITE(test_popins2){

    //SEQAN_CALL_TEST(test_ccdbg_opt);
    //SEQAN_CALL_TEST(test_ccdbg_build);
    //SEQAN_CALL_TEST(test_ccdbg_connected_components);
    //SEQAN_CALL_TEST(test_ccdbg_simpleBranching_singleThread);
    //SEQAN_CALL_TEST(test_ccdbg_simpleBubbles);
    //SEQAN_CALL_TEST(test_ccdbg_simpleColorTest);
    //SEQAN_CALL_TEST(test_ccdbg_rev_comp);

    SEQAN_CALL_TEST(setup_lecc_unittest);
    SEQAN_CALL_TEST(test_lecc_unittest);
}
SEQAN_END_TESTSUITE
