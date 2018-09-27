/*!
* \file    src/util.h
* \brief   Header for popins2 utility functions
*
*/
#ifndef POPINS2_UTIL_H_
#define POPINS2_UTIL_H_

#include <seqan/arg_parse.h>

#include <iostream>
#include <vector>
#include <algorithm>            // std::sort
#include <dirent.h>             // read folder



template <typename TType>
inline float median(const std::vector<TType> &vec){

    if(vec.size()==1) return vec[0];          // will happen a lot, catch first, no sort needed
    else if(vec.empty())  return 0;
    else{
        std::vector<uint16_t> v_copy(vec);

        size_t size = v_copy.size();

        std::sort(v_copy.begin(), v_copy.end());

        return size%2 == 0 ? (v_copy[size / 2 - 1] + v_copy[size / 2]) / 2 : v_copy[size / 2];
    }

}


inline void printTimeStatus(const char * message){
        // Get the current date and time.
        char timestamp[80];
        time_t now = time(0);
        struct tm tstruct;
        tstruct = *localtime(&now);
        strftime(timestamp, sizeof(timestamp), "[popins2 %Y-%m-%d %X] ", &tstruct);

        // Print time and message.
        std::cerr << timestamp << message << std::endl;
}


inline void printTimeStatus(std::ostringstream & message){
        std::string msg = message.str();
        printTimeStatus(seqan::toCString(msg));
}


/*!
* \fn       vector<string> getFilesFromDir(string &path)
* \brief    Function returnes a vector of all files in a given folder (path)
* \remark   Only works for UNIX so far.
* \return   vector<string>
*/
inline std::vector<std::string> getFilesFromDir(std::string &path){
    std::vector<std::string> sample_fastx_names;      // all file names in --indir with full path
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(path))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastx!="." && current_fastx!="..")
                sample_fastx_names.push_back(path+current_fastx);
        }
        closedir(dir);
    }
    else {
        // could not open directory
        std::perror ("");
    }

    return sample_fastx_names;
}


/*!
* \fn       void getFilesFromDir(std::vector<std::string> &v, std::string &path)
* \brief    Function fills a referenced vector with all files in a given folder (path).
* \remark   Only works for UNIX so far.
* \return   void
*/
inline void getFilesFromDir(std::vector<std::string> &v, std::string &path){

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(path))) != NULL){
        // print all the files and directories within directory
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;
            // exclude UNIX directory navigation links
            if (current_fastx!="." && current_fastx!="..")
                v.push_back(path+current_fastx);
        }
        closedir(dir);
    }
    else {
        // could not open directory
        std::perror ("");
    }
}




#endif /*POPINS2_UTIL_H_*/