/**
 *  @file   src/popins2_megamerge.h
 *  @brief  Creating a multi-k-iterative colored compacted de Bruijn Graph.
 *          In the past, the multi-k paradigm for de Bruijn Graphs (dBG) has been shown to be a successful
 *          method to transform the topology of a dBG for de novo sequence assembly. The multi-k paradigm
 *          aims to increase the average unitig length without the fragmentation of large k values.
 *          The megamerge module extends this idea for colored de Bruijn Graphs by maintaining the color
 *          awareness in each iteration.
 */

#ifndef POPINS2_MEGAMERGE_H_
#define POPINS2_MEGAMERGE_H_

#include <bifrost/ColoredCDBG.hpp>          // ColoredCDBG
#include "util.h"                           // getFastx, printTimeStatus
#include "argument_parsing.h"               // MegamergeOptions, parseCommandLine, printMegamergeOptions
#include "prettyprint.h"                    // REMOVE AT RELEASE, print


/**
 *          This function defines the boolean logic to decide if a color for unitig is significant.
 *  @param  b1, b2, b3 are boolean values to compare
 */
inline bool logicDecision(const bool b1, const bool b2, const bool b3){
    // true if at least 2 out of 3 are true
    return b1 ? (b2 || b3) : (b2 && b3);
}


/**
 *          This function extract the color bits of a kmer position of UnitigColorMap
 *  @param  ucm is a unitig of the graph
 *  @param  colorBits is a vector to store the color bits of the kmer
 *  @param  pos is the position of the kmer in the unitig (ucm)
 */
inline void getColorBitsOfPosition(const UnitigColorMap<void> &ucm, std::vector<bool> &colorBits, const size_t pos){
    UnitigColorMap<void> kmer;
    UnitigColors* colors;

    kmer   = ucm.getKmerMapping(pos);
    colors = kmer.getData()->getUnitigColors(kmer);

    UnitigColors::const_iterator cit = colors->begin(kmer);
    for (; cit != colors->end(); ++cit)
        colorBits[cit.getColorID()] = true;
}


/**
 *          This function defines how to subsample kmers from a unitig.
 *  @param  ucm is a unitig of the graph
 *  @param  color_indices is a list of samples the unitig should be added to
 */
inline void colorProbing(const UnitigColorMap<void> &ucm, std::vector<size_t> &color_indices, const size_t nb_colors){
    // the current approach is to ckeck the unitig colors at three positions: start, middle and end
    // NOTE: finding the end position usually requires strand awareness, but this is not important here
    size_t end = ucm.len - 1;
    size_t mid = floor(end / 2);

    std::vector<bool> startColorBits(nb_colors, false);
    std::vector<bool> middleColorBits(nb_colors, false);
    std::vector<bool> endColorBits(nb_colors, false);

    getColorBitsOfPosition(ucm, startColorBits, 0);
    getColorBitsOfPosition(ucm, middleColorBits, mid);
    getColorBitsOfPosition(ucm, endColorBits, end);

    for (size_t i = 0; i < nb_colors; ++i)
        if (logicDecision(startColorBits[i], middleColorBits[i], endColorBits[i]))
            color_indices.push_back(i);
}


/**
 *          Megamerge main routine.
 *  @return 0=success, 1=error
 */
int popins2_megamerge(int argc, char const *argv[]){
    // =====================
    // Argument parsing
    // =====================
    MegamergeOptions mmo;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(mmo, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "[popins2 megamerge] seqan::ArgumentParser::PARSE_ERROR" << endl;
        return 1;
    }
    printMegamergeOptions(mmo);

    std::ostringstream msg;
    msg << "[popins2 megamerge] Starting ..."; printTimeStatus(msg);

    std::vector<std::string> samples;
    getFastx(samples, toCString(mmo.samplePath));

    // manage a temporary directory
    if (strcmp(toCString(mmo.tempPath), "") == 0)   // if mmo.tempPath == ""
        (mkdir("megamerge_aux", 0750) == -1) ? cerr << "Error :  " << strerror(errno) << endl : cout << "Directory created." << "\n";

    int delta_k = 20;
    int k_max = 123;

    // =====================
    // Graph options
    // =====================
    CCDBG_Build_opt opt;
    opt.filename_ref_in   = samples;
    opt.deleteIsolated    = true;
    opt.clipTips          = true;
    opt.prefixFilenameOut = "ccdbg";
    opt.nb_threads        = 16;
    opt.outputGFA         = true;
    opt.verbose           = false;
    opt.k                 = 23;

    // =====================
    // Multi-k framework
    // =====================
    ColoredCDBG<> g(opt.k);

    msg << "[popins2 megamerge] Building dBG..."; printTimeStatus(msg);
    g.buildGraph(opt);

    msg << "[popins2 megamerge] Simplifying dBG..."; printTimeStatus(msg);
    g.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

    msg << "[popins2 megamerge] Annotating colors..."; printTimeStatus(msg);
    g.buildColors(opt);

    opt.k += delta_k;

    // exit point
    if (opt.k > k_max){
        msg << "[popins2 megamerge] Writing k_max GFA..."; printTimeStatus(msg);
        g.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
        return 0;
    }

    // get sample-unitig assignment
    const size_t nb_colors = g.getNbColors();
    std::unordered_map <std::string, std::vector<unsigned> > unitig_propagation_table;
    std::vector<size_t> color_indices;
    unsigned ucm_index = 0;

    for (auto &ucm : g){
        // don't process unitigs of size less than current k, they wouldn't survive the include of the next k iteration anyway
        if (ucm.size < (unsigned)opt.k){
            ++ucm_index;
            continue;
        }

        colorProbing(ucm, color_indices, nb_colors);

        for (size_t i = 0; i < color_indices.size(); ++i)
            unitig_propagation_table[g.getColorName(i)].push_back(ucm_index);

        ++ucm_index;
        color_indices.clear();
    }

    // DEBUG
    /*
    std::cout << "Table size: " << unitig_propagation_table.size() << '\n';
    for (auto cit = unitig_propagation_table.cbegin(); cit != unitig_propagation_table.cend(); ++cit){
        std::cout << cit->first << " : "; prettyprint::print(cit->second, 10);
    }
    */

    // write unitig file per sample
    for (auto sample = unitig_propagation_table.cbegin(); sample != unitig_propagation_table.cend(); ++sample){

        // open unitig file
        size_t lastSlashPos = sample->first.find_last_of("/");
        std::string fname   = sample->first.substr(lastSlashPos+1);
        size_t lastDotPos   = fname.find_last_of(".");
        std::string fname_  = fname.substr(0,lastDotPos);

        std::string temp_fname;
        if (strcmp(toCString(mmo.tempPath), "") == 0){      // if mmo.tempPath == ""
            temp_fname = getAbsoluteFileName("megamerge_aux", fname_) + ".unitigs.k" + std::to_string(opt.k) + ".fasta";
        }
        else{
            temp_fname = getAbsoluteFileName(toCString(mmo.tempPath), fname_) + ".unitigs.k" + std::to_string(opt.k) + ".fasta";
        }

        std::ofstream stream(temp_fname);
        if (!stream.good())
        {
            std::cerr << "ERROR: Could not open sample file \'" << temp_fname << "\' for writing." << std::endl;
            return 1;
        }

        // iterators pointing to the unitig ID list of a sample
        std::vector<unsigned>::const_iterator idx = sample->second.cbegin();
        std::vector<unsigned>::const_iterator idx_end = sample->second.cend();

        // loop through graph
        unsigned current_ucm_idx = 0;
        for (auto &ucm : g){

            if (*idx == current_ucm_idx){

                // write unitig to file
                stream << ">unitig_" << std::to_string(current_ucm_idx) << "\n";
                stream << ucm.referenceUnitigToString() << "\n";

                ++idx;

                if (idx == idx_end)
                    break;
            }

            current_ucm_idx += 1;
        }

        stream.close();
    }



    // =====================
    // EOF
    // =====================
    msg << "[popins2 megamerge] Done."; printTimeStatus(msg);

    return 0;
}



#endif /*POPINS2_MEGAMERGE_H_*/
