#include <iostream>
#include <bifrost/CompactedDBG.hpp>

int main()
{
    CDBG_Build_opt bf_options;

    bf_options.fastx_filename_in.push_back("/path/to/file/1.fq");
    bf_options.fastx_filename_in.push_back("/path/to/file/2.fq");
    bf_options.k = 31;
    bf_options.nb_unique_kmers = 24883500;
    bf_options.nb_non_unique_kmers = 8749400;
    bf_options.prefixFilenameOut = "myTestout";
    bf_options.nb_threads = 4;
    bf_options.clipTips = true;
    bf_options.deleteIsolated = true;

    // ---------- RUN CDBG ----------
    CompactedDBG<> cdbg(bf_options.k, bf_options.g);
    cout << "[PROGRESS] Building CDBG..." << endl;
    cdbg.build(bf_options);

    // ------- SIMPLIFY CDBG -------
    cout << "[PROGRESS] Simplifying CDBG..." << endl;
    cdbg.simplify(bf_options.deleteIsolated, bf_options.clipTips, bf_options.verbose);

    // --------- WRITE CDBG ---------
    cout << "[PROGRESS] Writing GFA..." << endl;
    cdbg.write(bf_options.prefixFilenameOut, bf_options.nb_threads, true, bf_options.verbose);

    cout << "The DBG has " << cdbg.size() << " unitigs.\n" << endl;
}
