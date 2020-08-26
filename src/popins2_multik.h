/**
 *  @file   src/popins2_multik.h
 *  @brief  Creating a multi-k-iterative colored compacted de Bruijn Graph.
 *          In the past, the multi-k paradigm for de Bruijn Graphs (dBG) has been shown to be a successful
 *          method to transform the topology of a dBG for de novo sequence assembly. The multi-k paradigm
 *          aims to increase the average unitig length without the fragmentation of large k values.
 *          The multik module extends this idea for colored de Bruijn Graphs by maintaining the color
 *          awareness in each iteration.
 */

#ifndef POPINS2_MULTIK_H_
#define POPINS2_MULTIK_H_

#include <bifrost/ColoredCDBG.hpp>          // ColoredCDBG
#include <seqan/seq_io.h>                   // getAbsolutePath, toCString, readRecords
#include "util.h"                           // getFastx, printTimeStatus, getAbsoluteFileName
#include "argument_parsing.h"               // MultikOptions, parseCommandLine, printMultikOptions

using namespace seqan;



/**
 *          Function to convert FASTQ to FASTA.
 * @details This is the seqAn2 way to do this. Converting this with the new C++20 concepts used in
 *          seqAn3 is so much more elegant.
 * @param   fastq_file is the name of the fastq file to convert to fasta, including its full path
 * @param   outpath is a (optional) string to define an output directory for the FASTA.
 *          If outpath is "" the FASTA will be written to the current working directory
 * @param   fasta_names is a vector reference to store the absolute FASTA filenames in
 * @return  bool; 1 if error, 0 else
 */
inline bool fastq2fastq(const std::string &fastq_file, const std::string &outpath, std::vector<std::string> &fasta_names){
    // read FASTQ
    CharString seqFileName = fastq_file;
    SeqFileIn seqFileIn;

    if (!open(seqFileIn, toCString(seqFileName))){
        std::cerr << "ERROR: Could not open FASTQ file to read from.\n";
        return 1;
    }

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    StringSet<CharString> quals;
    clear(quals);   // we don't need them

    try{
        readRecords(ids, seqs, quals, seqFileIn);
    }
    catch (Exception const & e){
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    close(seqFileIn);

    // get FASTQ filename without path and without file ending
    size_t lastSlashPos = fastq_file.find_last_of("/");
    std::string fname   = fastq_file.substr(lastSlashPos+1);
    size_t lastDotPos   = fname.find_last_of(".");
    std::string fname_  = fname.substr(0,lastDotPos);

    // create FASTA filename
    std::string ofile_name;
    if (strcmp(outpath.c_str(), "") == 0){
        ofile_name = fname_+".fasta";
    }
    else{
        ofile_name = getAbsoluteFileName(outpath, fname_+".fasta");
    }
    fasta_names.push_back(ofile_name);

    // write FASTA
    SeqFileOut seqFileOut(ofile_name.c_str());
    if (!open(seqFileOut, ofile_name.c_str())){
        std::cerr << "ERROR: Could not open FASTA file to write into.\n";
        return 1;
    }

    try{
        writeRecords(seqFileOut, ids, seqs);
    }
    catch (Exception const & e){
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }


    close(seqFileOut);

    return 0;
}


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
 *          Multik main routine.
 *  @return 0=success, 1=error
 */
int popins2_multik(int argc, char const *argv[]){
    // =====================
    // Argument parsing
    // =====================
    MultikOptions mko;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(mko, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK){
        if (res == seqan::ArgumentParser::PARSE_HELP)
            return 0;
        cerr << "[popins2 multik] seqan::ArgumentParser::PARSE_ERROR" << endl;
        return 1;
    }

    // manage a temporary directory
    if (strcmp(mko.tempPath.c_str(), "") != 0)      // if mko.tempPath != ""
        (mkdir(mko.tempPath.c_str(), 0750) == -1) ? cerr << "Error :  " << strerror(errno) << endl : cout << "[popins2 multik] Temp directory created.\n";

    int delta_k = mko.delta_k;
    int k_max   = mko.k_max;

    printMultikOptions(mko);

    std::ostringstream msg;
    msg << "[popins2 multik] Starting ..."; printTimeStatus(msg);

    // read original FASTQ samples
    std::vector<std::string> samples;
    getFastx(samples, mko.samplePath);

    // create a FASTA at tempDir for every input FASTQ
    std::vector<std::string> temp_fastas;
    bool fastq2fasta_failed = false;

    for (auto &sample : samples){
        bool ret = fastq2fastq(sample, mko.tempPath, temp_fastas);
        fastq2fasta_failed = fastq2fasta_failed || ret;
    }
    if (fastq2fasta_failed){
        std::cerr << "[Error] Initial FASTQ to FASTA conversion failed]" << '\n';
        return 1;
    }

    samples.clear();

    // =====================
    // Graph options
    // =====================
    CCDBG_Build_opt opt;
    opt.filename_seq_in   = temp_fastas;
    opt.deleteIsolated    = true;
    opt.clipTips          = true;
    opt.useMercyKmers     = false;
    opt.prefixFilenameOut = mko.prefixFilenameOut;
    opt.nb_threads        = 16;
    opt.outputGFA         = true;
    opt.verbose           = false;
    opt.k                 = mko.k_init;

    // =====================
    // Multi-k framework
    // =====================
    unsigned k_iter_counter = 0;
    while (opt.k <= k_max) {
        ++k_iter_counter;
        msg << "[popins2 multik] Multi-k iteration "+std::to_string(k_iter_counter)+" using k="+std::to_string(opt.k); printTimeStatus(msg);

        ColoredCDBG<> g(opt.k);

        msg << "[popins2 multik] Building dBG..."; printTimeStatus(msg);
        g.buildGraph(opt);

        msg << "[popins2 multik] Simplifying dBG..."; printTimeStatus(msg);
        g.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

        msg << "[popins2 multik] Annotating colors..."; printTimeStatus(msg);
        g.buildColors(opt);

        opt.k += delta_k;

        // exit point
        if (opt.k > k_max){
            msg << "[popins2 multik] Writing "+opt.prefixFilenameOut+".gfa of de Bruijn Graph built with k="+std::to_string(opt.k - delta_k)+"..."; printTimeStatus(msg);
            g.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
            break;
        }

        // get sample-unitig assignment
        msg << "[popins2 multik] Get sample-unitig assignment..."; printTimeStatus(msg);
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

        // write unitig file per sample
        msg << "[popins2 multik] Adding (k + delta_k)-unitigs to temp FASTAs..."; printTimeStatus(msg);
        for (auto sample = unitig_propagation_table.cbegin(); sample != unitig_propagation_table.cend(); ++sample){

            std::ofstream stream(sample->first, std::ofstream::out | std::ofstream::app);
            if (!stream.good())
            {
                std::cerr << "ERROR: Could not open sample file \'" << sample->first << "\' for writing." << std::endl;
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

    }

    // =====================
    // EOF
    // =====================
    msg << "[popins2 multik] Done."; printTimeStatus(msg);

    return 0;
}



#endif /*POPINS2_MULTIK_H_*/
