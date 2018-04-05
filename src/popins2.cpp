#include <iostream>
#include <ctime>

#include "argument_parsing.h"           /* seqAn argument parser */
#include "popins2_single.h"
#include "popins2_merge.h"

using namespace std;



// ==============================
// Function: main()
// ==============================
int main(int argc, char const *argv[]){



    //popins2_single(argc, argv);

    popins2_merge(argc, argv);



    return 0;
}



