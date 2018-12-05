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
#include <cerrno>


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


/*!
* \fn       void getFastx(std::vector<std::string> &v, std::string &path)
* \brief    Function fills a referenced vector with all FASTX (*.f[ast]?a|q[.gz]?) files in a given path.
* \remark   Only works for UNIX/POSIX so far.
* \return   bool; false if no fastx file was found
*/
inline bool getFastx(std::vector<std::string> &v, std::string &path, const bool verbose = false){
    bool ret = false;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(seqan::toCString(path))) != NULL){
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;

            if      (current_fastx=="." || current_fastx==".." || current_fastx=="contigs.fa")  continue;
            else if (current_fastx.substr(current_fastx.length()-1) == "/")                     continue;    // NOTE: edit here for recursive directory search
            else{
                // taken and adapted from: https://github.com/pmelsted/bifrost/blob/master/src/File_Parser.hpp
                const size_t last_point = current_fastx.find_last_of(".");
                std::string s_ext = current_fastx.substr(last_point + 1);

                if ((s_ext == "gz")){
                    s_ext = current_fastx.substr(current_fastx.find_last_of(".", last_point - 1) + 1);
                    if ((s_ext == "fasta.gz") || (s_ext == "fa.gz") || (s_ext == "fastq.gz") || (s_ext == "fq.gz")){
                        v.push_back(path+current_fastx);
                        ret = true;
                    }
                    else {
                        if (verbose){
                            std::cerr << "getFastx(): Compressed input file is not in FASTX (*.fasta.gz, *.fa.gz, *.fastq.gz, *.fq.gz) format." << std::endl;
                            std::cerr << "Skipping " << current_fastx << std::endl;
                        }
                    }
                }
                else if ((s_ext == "fasta") || (s_ext == "fa") || (s_ext == "fastq") || (s_ext == "fq")){
                    v.push_back(path+current_fastx);
                    ret = true;
                }
                else{
                    if (verbose){
                        std::cerr << "getFastx(): Input file is not in FASTX (*.f[ast]?a|q[.gz]?) format." << std::endl;
                        std::cerr << "Skipping " << current_fastx << std::endl;
                    }
                }
            }
        }
        closedir(dir);
    }
    else if (errno == ENOENT){
        // could not open directory
        std::perror ("No such file or directory.");
    }
    return ret;
}


inline bool file_exist (const std::string &name){
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}


#endif /*POPINS2_UTIL_H_*/
