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

    std::time_t start_time = std::time(0);

    int ret = 0;

    const char * prog_name = argv[0];
    if (argc < 2){
        printHelp(prog_name);
        return 1;
    }

    const char * command = argv[1];
    if (strcmp(command,"single") == 0) ret = popins2_single(argc, argv);
    else if (strcmp(command,"merge") == 0) ret = popins2_merge(argc, argv);
    else if (strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0){
        printHelp(prog_name);
        return 1;
    }
    else{
        std::cerr << "ERROR: Unknown command: " << command << std::endl;
        printHelp(prog_name);
        return 1;
    }

    if (ret == 1)
       return 1;

    if (ret == 0){
        std::ostringstream msg;
        msg << "popins2 " << command << " finished in " << (std::time(0) - start_time) << " seconds.";
        printTimeStatus(msg);
    }

    return 0;

}



