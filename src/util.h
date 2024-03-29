/*!
* \file    src/util.h
* \brief   Header for popins2 utility functions
*
*/
#ifndef POPINS2_UTIL_H_
#define POPINS2_UTIL_H_

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <vector>
#include <algorithm>            // std::sort
#include <dirent.h>             // read folder
#include <cerrno>

using namespace seqan;



inline std::string forcePathEndingSlash(std::string path){
    if (path.substr(path.length()-1) != "/")
        path += "/";
    return path;
}


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
        message.str("");
}


inline bool file_exist (const std::string &name){
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
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
* \fn       void getFastx(std::vector<std::string> &v, std::string path, const bool verbose = false)
* \brief    Function fills a referenced vector with all FASTX (*.f[ast]?a|q[.gz]?) files in a given path.
* \remark   Only works for UNIX/POSIX so far.
* \return   bool; false if no fastx file was found
*/
inline bool getFastx(std::vector<std::string> &v, std::string path, const bool verbose = false){
    bool ret = false;
    DIR *dir;
    struct dirent *ent;

    if (path.back() != '/')
        path += "/";

    if ((dir = opendir(seqan::toCString(path))) != NULL){
        while ((ent = readdir(dir)) != NULL){
            std::string current_fastx = ent->d_name;

            if      (current_fastx=="." || current_fastx=="..")                                 continue;
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
                        std::cerr << "getFastx(): Input file is not in FASTX (*.f[ast]?[a|q][.gz]?) format." << std::endl;
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


/*!
* \fn       listFiles(vector<string> &v, const string &path, const string &contigsFileName)
* \brief    Lists all files <prefix> --> <sample> --> <filename>
* \remark   Only works for UNIX/POSIX so far.
* \return   void
*/
inline void listFiles(std::vector<std::string> &v, const std::string &path, const std::string &contigsFileName){

    DIR *dir = opendir(path.c_str());

    struct dirent *entry = readdir(dir);
    while (entry != NULL){

        if (entry->d_type == DT_DIR){

            std::string sampleID = entry->d_name;

            if (sampleID == "." || sampleID == ".."){
                entry = readdir(dir);
                continue;
            }

            std::string fileWithPath;
            if (path.back() != '/')
                fileWithPath = path+"/"+sampleID+"/"+contigsFileName;
            else
                fileWithPath = path+sampleID+"/"+contigsFileName;


            if (file_exist(fileWithPath))
                v.push_back(fileWithPath);
            else
                std::cerr << "WARNING: \'" << fileWithPath << "\' does not exist." << std::endl;
        }
        entry = readdir(dir);
    }

    closedir(dir);
}


inline std::string getAbsoluteFileName(const std::string &path, const std::string &file)
{
    std::string abspath(path);

    if (abspath.substr(abspath.length()-1) != "/")
        abspath += "/";

    abspath += file;

    return abspath;
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// copied from popins-1.0.1
// --------------------------------------------------------------------------------------------------------------------------------------------

// ==========================================================================
// Struct SampleInfo
// ==========================================================================

struct SampleInfo
{
    CharString sample_id;
    CharString bam_file;
    double avg_cov;
    unsigned read_len;
    CharString adapter_type;

    SampleInfo() {}
};

// ==========================================================================
// Function initSampleInfo()
// ==========================================================================

SampleInfo
initSampleInfo(CharString & filename, CharString sample_id, CharString & adapter_type)
{
    SampleInfo info;
    info.bam_file = filename;
    info.sample_id = sample_id;

    BamFileIn bamFile(toCString(filename));
    BamHeader header;
    readHeader(header, bamFile);
    BamAlignmentRecord record;
    readRecord(record, bamFile);
    info.read_len = length(record.seq);

    info.adapter_type = adapter_type;

    return info;
}

// ==========================================================================
// Function readSampleInfo()
// ==========================================================================

bool
readSampleInfo(SampleInfo & info, CharString & filename)
{
    std::ifstream stream(toCString(filename));
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open sample info file \'" << filename << "\' for reading." << std::endl;
        return 1;
    }

    std::string field, value;
    while (stream >> field >> value)
    {
        if (field.compare("SAMPLE_ID") == 0)
            info.sample_id = value;
        else if (field.compare("BAM_FILE") == 0)
            info.bam_file = value;
        else if (field.compare("AVG_COV") == 0)
            info.avg_cov = lexicalCast<double>(value);
        else if (field.compare("READ_LEN") == 0)
            lexicalCast<unsigned>(info.read_len, value);
        else if (field.compare("ADAPTER_TYPE") == 0)
            info.adapter_type = value;
        else
            std::cerr << "WARNING: Ignoring field \'" << field << "\' in sample info file \'" << filename << "\'." << std::endl;
    }

    return 0;
}

// ==========================================================================
// Function write()
// ==========================================================================

bool
writeSampleInfo(SampleInfo & info, CharString & filename)
{
    // open the fileDate
    std::ofstream stream(toCString(filename));
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open sample info file \'" << filename << "\' for writing." << std::endl;
        return 1;
    }

    stream << "SAMPLE_ID" << "\t" << info.sample_id << "\n";
    stream << "BAM_FILE" << "\t" << info.bam_file << "\n";
    stream << "AVG_COV" << "\t" << info.avg_cov << "\n";
    stream << "READ_LEN" << "\t" << info.read_len << "\n";
    stream << "ADAPTER_TYPE" << "\t" << info.adapter_type << "\n";

    stream.close();
    return 0;
}

// ==========================================================================

void printStatus(const char * message)
{
        // Get the current date and time.
        char timestamp[80];
        time_t now = time(0);
        struct tm tstruct;
        tstruct = *localtime(&now);
        strftime(timestamp, sizeof(timestamp), "[popins2 %Y-%m-%d %X] ", &tstruct);

        // Print time and message.
        std::cerr << timestamp << message << std::endl;
}

void printStatus(std::ostringstream & message)
{
        std::string msg = message.str();
        printStatus(toCString(msg));
}

// ==========================================================================

// Returns true if file exists, otherwise false.
inline bool exists(CharString const & filename)
{
    struct stat buffer;
    return (stat(toCString(filename), &buffer) == 0);
}

// ==========================================================================

// Lists all files <prefix>/*/<filename>.
String<Pair<CharString> >
listFiles(CharString & prefix, CharString & filename)
{
   String<Pair<CharString> > paths;

    DIR *dir = opendir(toCString(prefix));

    struct dirent *entry = readdir(dir);
    while (entry != NULL)
    {
        if (entry->d_type == DT_DIR)
        {
           CharString sampleID = entry->d_name;
           if (sampleID == "." || sampleID == "..")
           {
               entry = readdir(dir);
               continue;
           }
           std::stringstream path;
           path << prefix << "/" << sampleID << "/" << filename;
           CharString pathStr = path.str();
           if (exists(pathStr))
                appendValue(paths, Pair<CharString>(sampleID, pathStr));
           else
              std::cerr << "WARNING: \'" << pathStr << "\' does not exist." << std::endl;
        }
        entry = readdir(dir);
    }

    closedir(dir);

    return paths;
}

// ==========================================================================

// List all directories <prefix>/*/ and return only the basename.
String<CharString>
listSubdirectories(CharString & prefix)
{
   String<CharString> subdirs;

   DIR *dir = opendir(toCString(prefix));

   struct dirent *entry = readdir(dir);
   while (entry != NULL)
   {
      if (entry->d_type == DT_DIR)
      {
          CharString sampleID = entry->d_name;
          if (sampleID == "." || sampleID == "..")
          {
              entry = readdir(dir);
              continue;
          }
         appendValue(subdirs, sampleID);
      }
      entry = readdir(dir);
   }

   closedir(dir);

   return subdirs;
}

// ==========================================================================

inline CharString
getFileName(CharString const & path, CharString const & name)
{
    CharString filename = path;
    filename += "/";
    filename += name;
    return filename;
}

// ==========================================================================

inline void
removeFile(CharString const & path, const char * filename)
{
    CharString file = path;
    file += "/";
    file += filename;
    remove(toCString(file));
}

// ==========================================================================
// Function checkFileEnding()
// ==========================================================================

bool
checkFileEnding(CharString & filename, std::string ending)
{
    std::string name = toCString(filename);
    size_t dotPos = name.find_last_of('.');

    if (dotPos == std::string::npos)
        return false;

    return name.substr(dotPos + 1, 3) == ending;
}

// ==========================================================================
// Function readFileNames()
// ==========================================================================

bool
readFileNames(String<CharString> & files, CharString & filenameFile)
{
    if (filenameFile == "") return 0;

    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    std::string file;
    while (stream >> file)
        appendValue(files, CharString(file));

    return 0;
}

// ==========================================================================

template <typename TValue>
bool readFileNames(String<CharString> & files, String<TValue> & values)
{
    if (length(files) > 1) return 0;
    std::cerr << "ReadFileNames " << length(files) << " " << files << std::endl;
    // Open input file
    CharString filenameFile = files[0];
    std::cerr << filenameFile << " " << files << std::endl;
    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    clear(files);
    clear(values);

    std::string file;
    TValue val;
    while (stream >> file >> val)
    {
        appendValue(files, file);
        appendValue(values, val);
    }

    return 0;
}

// ==========================================================================

bool
parseInterval(Triple<CharString, unsigned, unsigned> & out, CharString & in)
{
    Iterator<CharString, Rooted>::Type it = begin(in, Rooted());

    unsigned colonPos = 0;
    while (it != end(in))
    {
        if (*it == ':')
        {
            colonPos = position(it);
            break;
        }
        ++it;
    }

    if (colonPos == 0)
    {
        out.i1 = in;
        out.i2 = 0;
        out.i3 = maxValue<unsigned>();
        return 0;
    }

    unsigned dashPos = 0;
    while (it != end(in))
    {
        if (*it == '-')
        {
            dashPos = position(it);
            break;
        }
        ++it;
    }

    if (dashPos == 0)
    {
        std::cerr << "ERROR: Interval is not in format CHR:BEG-END." << std::endl;
        return 1;
    }

    out.i1 = prefix(in, colonPos);
    out.i2 = lexicalCast<unsigned>(infix(in, colonPos + 1, dashPos));
    out.i3 = lexicalCast<unsigned>(suffix(in, dashPos + 1));

    return 0;
}

bool
readChromosomes(std::set<CharString> & chromosomes, CharString & referenceFile)
{
    // Load or build and save the FASTA index.
    FaiIndex faiIndex;
    if (!open(faiIndex, toCString(referenceFile)))
    {
        if (!build(faiIndex, toCString(referenceFile)))
        {
            std::cerr << "ERROR: FASTA index could not be loaded or built.\n";
            return 1;
        }
        if (!save(faiIndex))    // Name is stored from when reading.
        {
            std::cerr << "WARNING: FASTA index could not be written to disk.\n";
        }
    }

    for (unsigned i = 0; i < length(faiIndex.seqNameStore); ++i)
        chromosomes.insert(faiIndex.seqNameStore[i]);

    return 0;
}




#endif /*POPINS2_UTIL_H_*/
