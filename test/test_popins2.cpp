#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/seq_io.h>

#include <bifrost/ColoredCDBG.hpp>
#include <../src/ColoredDeBruijnGraph.h>
#include <../src/LECC_Finder.h>

typedef std::unordered_map<Kmer, bool, KmerHash> border_map_t;

typedef uint8_t direction_t;

typedef std::unordered_map<uint64_t, Kmer> jump_map_t;

const direction_t VISIT_SUCCESSOR   = 0x0;
const direction_t VISIT_PREDECESSOR = 0x1;

using namespace seqan;

class LECC_Finder_Tester{

    LECC_Finder* f_;

public:

    LECC_Finder_Tester(LECC_Finder* f) : f_(f) {}

    void test_get_borders(border_map_t &m, const unsigned nb_leccs) const{
        f_->get_borders(m, nb_leccs);
    }

    void test_check_accessibility(border_map_t &m, const Kmer &kmer, const unsigned nb_leccs) const{
        f_->check_accessibility(m, kmer, nb_leccs);
    }

};

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


template <typename TKey, typename TValue>
inline void print(std::unordered_map<TKey, TValue> &m){
    // stay safe with constant iterators
    typename std::unordered_map<TKey,TValue>::const_iterator it;
    for (it = m.cbegin(); it != m.cend(); ++it){
        // output format
        std::cout << "[" << it->first << "\t:\t" << it->second << "]" << std::endl;
    }
}


inline void print_borders(const border_map_t &border_kmers, ExtendedCCDBG &xg){

    // the loop double-checks the right Kmers with corresponding colors
    for (border_map_t::const_iterator cit = border_kmers.cbegin(); cit != border_kmers.cend(); ++cit){

        UnitigColorMap<UnitigExtension> ucm = xg.find(cit->first, true);

        const UnitigColors* kmer_colors = ucm.getData()->getUnitigColors(ucm);

        std::vector<size_t> color_ids;

        UnitigColors::const_iterator cit_ = kmer_colors->begin(ucm);
        for (; cit_ != kmer_colors->end(); ++cit_)
            color_ids.push_back(cit_.getColorID());

        std::string kmer_seq = cit->first.toString();
        bool accessibility_bit = cit->second;
        //std::string ucm_seq  =  ucm.mappedSequenceToString();     // should be the same output as the line before
        std::cout << "BORDER [" << kmer_seq << " : " << accessibility_bit << "] OF COLORS "; print(color_ids);
    }

}


inline void print_jump_map(const jump_map_t &m){
    for (jump_map_t::const_iterator cit = m.cbegin(); cit != m.cend(); ++cit)
        cout << cit->first << " -> " << cit->second.toString() << endl;
    cout << endl;
}


inline void print_unitig_ends(ExtendedCCDBG &g){
    for (auto &ucm : g){

        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* data = da->getData(ucm);

        size_t len = ucm.len;   // I assume this gets me the past-last-kmer index
        std::cout << data->getID() << " consists of " << len << "kmers." << std::endl;

        UnitigColorMap<UnitigExtension> k_first = ucm.getKmerMapping(0);    // seems to be 31 bp long, expected 63 bp
        UnitigColorMap<UnitigExtension> k_last  = ucm.getKmerMapping(len-1);    // seems to be 31 bp long, expected 63 bp

        const UnitigColors* k_first_colors = k_first.getData()->getUnitigColors(k_first);
        const UnitigColors* k_last_colors  =  k_last.getData()->getUnitigColors(k_last);

        std::vector<size_t> k_first_color_ids;
        std::vector<size_t> k_last_color_ids;

        // get color IDs of unitig's first kmer
        UnitigColors::const_iterator cit = k_first_colors->begin(k_first);
        for (; cit != k_first_colors->end(); ++cit){
            k_first_color_ids.push_back(cit.getColorID());
        }

        // get color IDs of unitig's last kmer
        cit = k_last_colors->begin(k_last);
        for (; cit != k_last_colors->end(); ++cit){
            k_last_color_ids.push_back(cit.getColorID());
        }

        // PRINT
        std::string first_seq = k_first.mappedSequenceToString();
        std::string last_seq  =  k_last.mappedSequenceToString();
        std::cout << data->getID() << " START COLORS (" << first_seq << "): " ; print(k_first_color_ids);
        std::cout << data->getID() << " END COLORS ("   << last_seq  << "): " ; print(k_last_color_ids);
    }
}


inline void init_ids_as_in_schematic_lecc_unittest(ExtendedCCDBG &xg){
    const char* u1_head = "CGTCAGCCGTCACCGAGCCTTATCTATGTCCCGACCAACAGAATAAACAGCAGCTACGCTTGA";
    Kmer u1_kmer = Kmer(u1_head);
    UnitigColorMap<UnitigExtension> u1 = xg.find(u1_kmer, true);
    DataAccessor<UnitigExtension>* da1 = u1.getData();
    UnitigExtension* data1 = da1->getData(u1);
    data1->setID(1u);

    const char* u2_head = "ATATATACATTTATAAATTTATATAAATATATATATTTATATATATATTTATCTTTATATTTG";
    Kmer u2_kmer = Kmer(u2_head);
    UnitigColorMap<UnitigExtension> u2 = xg.find(u2_kmer, true);
    DataAccessor<UnitigExtension>* da2 = u2.getData();
    UnitigExtension* data2 = da2->getData(u2);
    data2->setID(2u);

    const char* u3_head = "ACTCCGTGAATTGATCTGTTTGTCGAATACCCATGGTGTAAAGTGGAGGGTCCTACCCGCAAT";
    Kmer u3_kmer = Kmer(u3_head);
    UnitigColorMap<UnitigExtension> u3 = xg.find(u3_kmer, true);
    DataAccessor<UnitigExtension>* da3 = u3.getData();
    UnitigExtension* data3 = da3->getData(u3);
    data3->setID(3u);

    const char* u4_head = "ATATATACATTTATAAATTTATATAAATATATATATTTATATATATATTTATCTTTATATTTC";
    Kmer u4_kmer = Kmer(u4_head);
    UnitigColorMap<UnitigExtension> u4 = xg.find(u4_kmer, true);
    DataAccessor<UnitigExtension>* da4 = u4.getData();
    UnitigExtension* data4 = da4->getData(u4);
    data4->setID(4u);

    const char* u5_head = "TATATATACATTTATAAATTTATATAAATATATATATTTATATATATATTTATCTTTATATTT";
    Kmer u5_kmer = Kmer(u5_head);
    UnitigColorMap<UnitigExtension> u5 = xg.find(u5_kmer, true);
    DataAccessor<UnitigExtension>* da5 = u5.getData();
    UnitigExtension* data5 = da5->getData(u5);
    data5->setID(5u);
}


inline void simple_lecc_unittest_truthset(/*TODO*/){

    std::string t1 = "ACTCCGTGAATTGATCTGTTTGTCGAATACCCATGGTGTAAAGTGGAGGGTCCTACCCGCAATAAAGTACCCTATCGCTAAATACCTATTGACGTGTTTTGACGGAATTTAAACTACAACGGCGGTGACATTAATGAGCTGAAGCCCGACCGCCCGCTGGGCGCGCCTCCTTGGGCAAGTTTTTTTTTTTTTTTTTTTTTGCTTTTAGGCCACATTCGGGTGGCGTCAACAGAAACGGTATATATACATTTATAAATTTATATAAATATATATATTTATATATATATTTATCTTTATATTTCCGGAGTTCAAATGAAAGATGCACAGCATTCACTCTATTTTTTTTTTTTTTTTTTTTTGCAACTGAATGCACCGATATAGTGCGTACCGCGTTGTATTCTTCTTGACAGTAAGTAATCCGTTCCATCGGTCAAGCTAGACCACTCGAACGTTTGAAGTGCTTTTGGTACCCTTTCACGGAGTTCACAAGATGGGCGCCTGGTTCGTTTGTGTTTGGGAAGGGCGGTGTCCATAGGCCA";

    std::string t2 = "CGTCAGCCGTCACCGAGCCTTATCTATGTCCCGACCAACAGAATAAACAGCAGCTACGCTTGAGGAAGCACAAGTAGCGTCCCAGTGGTTACTCCAATTACTGCCTGTATTGACTAAAGCAGCAGTGGTCGCCTCTACCCAAACAGGACATCCATAACTGTCTGAGCACTGGGTGTCACTGGGGGGGGGGGGGGGGGGGGTCACGACCAGAAACTAAGTAATAAAGGCGTTCGGCCTATATATATACATTTATAAATTTATATAAATATATATATTTATATATATATTTATCTTTATATTTGCACACCTTCCTTCGGAACGACCGTCAGATTCACCTGCGGGGGGGGGGGGGGGGGGGGAAGGAAGCCGAGAACTCAAACCCACGCACCCTATGAGTAGTCCTACACGTCCATCGAGACCTCTCGTACTTGATGACACACACCCATATGCCAAATTCTTATCTGTACCTCTTATTCAAGAAGCGTGGCCCGACCTCACATCATATACACTGAGCATGGTGGTCGATATTACAGAACATG";

}


// --------------------
// | SIMPLE LECC TEST |
// --------------------
CCDBG_Build_opt opt_lecc_unittest;

SEQAN_DEFINE_TEST(setup_simple_lecc_unittest){

    std::vector<std::string> i_files;

    i_files.push_back("./data/lecc_unittest/sample1.fa");
    i_files.push_back("./data/lecc_unittest/sample2.fa");
    i_files.push_back("./data/lecc_unittest/sample3.fa");
    i_files.push_back("./data/lecc_unittest/sample4.fa");

    opt_lecc_unittest.filename_ref_in = i_files;
    opt_lecc_unittest.deleteIsolated = true;
    opt_lecc_unittest.clipTips = true;
    opt_lecc_unittest.prefixFilenameOut = "lecc_unittest";
    opt_lecc_unittest.nb_threads = 1;
    opt_lecc_unittest.outputGFA = true;
    opt_lecc_unittest.verbose = false;
    opt_lecc_unittest.k = 63;
}


SEQAN_DEFINE_TEST(simple_lecc_unittest){

    ExtendedCCDBG xg(opt_lecc_unittest.k, opt_lecc_unittest.g);

    /* DEBUG */ std::cout << "MAX_KMER_SIZE=" << MAX_KMER_SIZE << std::endl;
    /* DEBUG */ std::cout << "ExtendedCCDBG::getK()=" << xg.getK() << std::endl << std::endl;

    SEQAN_ASSERT_EQ(xg.buildGraph(opt_lecc_unittest),true);

    SEQAN_ASSERT_EQ(xg.simplify(opt_lecc_unittest.deleteIsolated, opt_lecc_unittest.clipTips, opt_lecc_unittest.verbose),true);

    SEQAN_ASSERT_EQ(xg.buildColors(opt_lecc_unittest),true);

    xg.init_ids();

    init_ids_as_in_schematic_lecc_unittest(xg);

    //std::cout << "---------- ALL UNITIG ENDS ----------" << std::endl;
    //print_unitig_ends(xg); std::cout << std::endl;

    xg.init_entropy();

    LECC_Finder F(&xg, 0.7f);

    unsigned nb_leccs = F.annotate();

    SEQAN_ASSERT_EQ(nb_leccs, 1u);      // simple_lecc_unittest has only one LECC

    LECC_Finder_Tester T(&F);

    border_map_t border_kmers;                 // storage for the borders per LECC
    T.test_get_borders(border_kmers, 1u);

    SEQAN_ASSERT_EQ(border_kmers.size(), 4u);   // simple_lecc_unittest has 4 border kmers

    std::unordered_map<uint64_t, unsigned> accessibility_truth_set({
        {17419696912337850218u, 2},
        {5617511096006481569u, 2},
        {16311619577470620275u, 2},
        {6202449885273403341u, 2}
    });

    // for every border kmer check if the correct amount of partners was accessed
    for (auto &border : border_kmers){

        T.test_check_accessibility(border_kmers, border.first, 1u);

        unsigned counter = 0;
        for (border_map_t::const_iterator cit = border_kmers.cbegin(); cit != border_kmers.cend(); ++cit)
            if (true == cit->second)
                ++counter;
        SEQAN_ASSERT_EQ(counter, 2u);     // in simple_lecc_unittest each border should be able to access 2 partners

        size_t ret = accessibility_truth_set.erase(border.first.hash());
        SEQAN_ASSERT_EQ(ret, 1u);   // erase() should have removed one element from accessibility_truth_set

        // reset accessibility bits
        for (auto &border : border_kmers) border.second = false;
    }

    SEQAN_ASSERT_EQ(accessibility_truth_set.size(), 0u);    // if all lecc borders were processed and the correct hash value was used, then accessibility_truth_set should be empty at this point

    jump_map_t jump_map;
    bool ret0 = F.find_jumps(jump_map, 1u);
    SEQAN_ASSERT_EQ(ret0, true);    // check for seccessful execution

    std::cout << "---------- ALL JUMP PAIRS ----------" << std::endl;
    print_jump_map(jump_map);

    SEQAN_ASSERT_EQ(xg.write(opt_lecc_unittest.prefixFilenameOut, opt_lecc_unittest.nb_threads, opt_lecc_unittest.verbose), true);
}


// -----------------------
// | 5 SIMULATED SIMPLES |
// -----------------------
CCDBG_Build_opt opt_5simu_test;

SEQAN_DEFINE_TEST(setup_5simu_test){

    std::vector<std::string> i_files;

    i_files.push_back("./data/5_simu_samples/S0001_assembly_k141.contigs.fa");
    i_files.push_back("./data/5_simu_samples/S0002_assembly_k141.contigs.fa");
    i_files.push_back("./data/5_simu_samples/S0003_assembly_k141.contigs.fa");
    i_files.push_back("./data/5_simu_samples/S0004_assembly_k141.contigs.fa");
    i_files.push_back("./data/5_simu_samples/S0005_assembly_k141.contigs.fa");

    opt_5simu_test.filename_ref_in = i_files;
    opt_5simu_test.deleteIsolated = true;
    opt_5simu_test.clipTips = true;
    opt_5simu_test.prefixFilenameOut = "5simu_test";
    opt_5simu_test.nb_threads = 1;
    opt_5simu_test.outputGFA = true;
    opt_5simu_test.verbose = false;
    opt_5simu_test.k = 31;
}


SEQAN_DEFINE_TEST(call_5simu_test){

    ExtendedCCDBG xg(opt_5simu_test.k, opt_5simu_test.g);

    /* DEBUG */ std::cout << "MAX_KMER_SIZE=" << MAX_KMER_SIZE << std::endl;
    /* DEBUG */ std::cout << "ExtendedCCDBG::getK()=" << xg.getK() << std::endl << std::endl;

    SEQAN_ASSERT_EQ(xg.buildGraph(opt_5simu_test),true);

    SEQAN_ASSERT_EQ(xg.simplify(opt_5simu_test.deleteIsolated, opt_5simu_test.clipTips, opt_5simu_test.verbose),true);

    SEQAN_ASSERT_EQ(xg.buildColors(opt_5simu_test),true);

    xg.init_ids();

    //std::cout << "---------- ALL UNITIG ENDS ----------" << std::endl;
    //print_unitig_ends(xg); std::cout << std::endl;

    xg.init_entropy();

    LECC_Finder F(&xg, 0.7f);

    unsigned nb_leccs = F.annotate();

    SEQAN_ASSERT_EQ(nb_leccs, 9u);

    LECC_Finder_Tester T(&F);

    std::unordered_map<std::string, unsigned> nb_border_kmers_per_lecc_truthset({       // one entry is a representing kmer for a border set
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", 23},
        {"AATATATAAATAAATATTATATATTATATAT", 1},
        {"TAGATGGGTGGATGGGTAGATGGGTGGATGA", 2},
        {"ACACACACACACACACACACACACACTTCTC", 16},
        {"ACAGTCTCAGACACAGACACACACCACACAC", 2},
        {"AGAGCGAGACTCCATCTCAAAAAAAAAAAAA", 4},
        {"ACATGTGCATACATGTATATACACATATGTG", 2},
        {"GTAATAATAATAATAATAATAATAATAATAA", 6},
        {"ATCTATCTATCTATCTATCTATCTATCTATC", 4}
    });

    std::unordered_map<std::string, unsigned> border_kmers_accessibility_truthset({
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGG", 8},
        {"ATTTTTTTTTTTTTGAGATGGAGTTTTGCTC", 0},
        {"CTTTTTTTTTTTTTGAGATGGAGTTTTGCTC", 0},
        {"TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 4},
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGA", 8},
        {"TCTTTTTTTTTTTTTTTTTTTTTTTTGAGAA", 2},
        {"GTCTCAAAAAAAAAAAAAAAAAAAAAAAGAC", 5},
        {"CTCTTTTTTTTTTTTTTTTTTTTTTTGAGAC", 5},
        {"AAAAAAAAAAAAAAAAAAAAAAAAGAAAAAA", 3},
        {"GACTCCATCTCAAAAAAAAAAAAAAAAAAAA", 4},
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT", 8},
        {"GTCTTTTTTTTTTTTTTTTTTTTTTTTGAGA", 2},
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", 8},
        {"CTTTTTTTTTTTTTTTTTTTTTTTTTGAGAC", 5},
        {"ACTCCGTCTCAAAAAAAAAAAAAAAAAAAAG", 2},
        {"AACTCTGTCTCAAAAAAAAAAAAAAAAAAAA", 7},
        {"GAGACTCCGTCTCAAAAAAAAAAAAAAAAAA", 11},
        {"GGAGCAAAACTCCATCTCAAAAAAAAAAAAA", 4},
        {"AAGACTCCGTCTCAAAAAAAAAAAAAAAAAA", 11},
        {"ATTTTTTTTTTTTTTTTTGAGACAGAGTCTC", 1},
        {"GACAGAGTGAGACTCTGTCTCAAAAAAAAAA", 8},
        {"AGAGACTCTGTCTCAAAAAAAAAAAAAAAAA", 7},
        {"AAAAAAAAAAAAAAAAAAAAAAAAGAAAAAT", 3},
        {"AATATATAAATAAATATTATATATTATATAT", 0},
        {"ATCTACCCATCCACCCATCTACCCATCCACA", 1},
        {"TAGATGGGTGGATGGGTAGATGGGTGGATGA", 1},
        {"GAACACACACACACACACACACACACACACA", 8},
        {"TCTCACACACACACACACACACACACACACA", 8},
        {"CCTCACACACACACACACACACACACACACA", 8},
        {"AAACACACACACACACACACACACACACACA", 8},
        {"CACACACACACACACACACACACACACACTA", 8},
        {"GACACACACACACACACACACACACACACAC", 8},
        {"ACACACACACACACACACACACACACACACG", 8},
        {"ACACACACACACACACACACACACACACACC", 8},
        {"ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTG", 8},
        {"ACACACACACACACACACACACACACTTCTG", 8},
        {"ACACACACACACACACACACACACACTTCTC", 8},
        {"CATCACACACACACACACACACACACACACA", 8},
        {"AATCACACACACACACACACACACACACACA", 8},
        {"GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTA", 8},
        {"CACACACACACACACACACACACACACACAG", 8},
        {"CACACACACACACACACACACACACACACAA", 8},
        {"AGTGTGTGGTGTGTGTCTGTGTCTGAGACTG", 1},
        {"ACAGTCTCAGACACAGACACACACCACACAC", 1},
        {"GACTCCATCTCAAAAAAAAAAAAAAGAAAAC", 2},
        {"GACTCCATCTCAAAAAAAAAAAAAAGAAAAA", 2},
        {"AGAGCGAGACTCCATCTCAAAAAAAAAAAAA", 2},
        {"TAGACTCCATCTCAAAAAAAAAAAAAAGAAA", 2},
        {"ACATGTATACACATGTACATATATACATATA", 1},
        {"ACATGTGCATACATGTATATACACATATGTG", 1},
        {"AATAATAATAATAATAATAATAATAATAAAT", 3},
        {"CAAATAATAATAATAATAATAATAATAATAA", 3},
        {"AATAATAATAATAATAATAATAATAATAAAG", 3},
        {"GTAATAATAATAATAATAATAATAATAATAA", 3},
        {"TAAATAATAATAATAATAATAATAATAATAA", 3},
        {"AATAATAATAATAATAATAATAATAATAAAA", 3},
        {"ATCTATCTATCTATCTATCTATCTATCTATC", 2},
        {"ATATCTATCTATCTATCTATCTATCTATCTA", 2},
        {"TAGATAGATAGATAGATAGATAGATAGATAA", 2},
        {"ATCTATCTATCTATCTATCTATCTATCTATA", 2}
    });

    for (size_t i = 1; i <= nb_leccs; ++i) {
        border_map_t border_kmers;                 // storage for the borders per LECC
        T.test_get_borders(border_kmers, i);
        cout << "border_kmers size " << border_kmers.size() << endl; cout << endl;

        // TEST correct amount of border kmers per lecc
        for (auto &t : nb_border_kmers_per_lecc_truthset){
            const char* s = t.first.c_str();
            Kmer sk(s);
            if (border_kmers.find(sk) != border_kmers.end()){
                SEQAN_ASSERT_EQ(border_kmers.size(), t.second);
                break;
            }
        }

        // TEST correct amount of accessible partners per border kmer
        for (auto &border : border_kmers){

            T.test_check_accessibility(border_kmers, border.first, i);

            unsigned counter = 0;
            for (border_map_t::const_iterator cit = border_kmers.cbegin(); cit != border_kmers.cend(); ++cit)
                if (true == cit->second)
                    ++counter;

            SEQAN_ASSERT_EQ(border_kmers_accessibility_truthset[border.first.toString()], counter);

            // reset accessibility bits
            for (auto &border : border_kmers) border.second = false;
        }
    }

    jump_map_t jump_map;
    bool ret0 = F.find_jumps(jump_map, nb_leccs);
    SEQAN_ASSERT_EQ(ret0, true);    // check for seccessful execution

    std::cout << "---------- ALL JUMP PAIRS ----------" << std::endl;
    print_jump_map(jump_map); cout << endl;

    SEQAN_ASSERT_EQ(xg.write(opt_5simu_test.prefixFilenameOut, opt_5simu_test.nb_threads, opt_5simu_test.verbose), true);
}


// --------------
// | CALL TESTS |
// --------------
SEQAN_BEGIN_TESTSUITE(test_popins2){

    SEQAN_CALL_TEST(setup_simple_lecc_unittest);

    SEQAN_CALL_TEST(simple_lecc_unittest);

    SEQAN_CALL_TEST(setup_5simu_test);

    SEQAN_CALL_TEST(call_5simu_test);
}


SEQAN_END_TESTSUITE
