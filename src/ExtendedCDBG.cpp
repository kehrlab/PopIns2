#include "ExtendedCDBG.h"
#include "CDBG_Data_extension.h"
#include <seqan/misc/union_find.h>



ExtendedCDBG::ExtendedCDBG(int kmer_length, int minimizer_length): CompactedDBG< UnitigExtension >(kmer_length, minimizer_length){
}

