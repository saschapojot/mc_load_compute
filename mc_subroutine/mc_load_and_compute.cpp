#include "mc_load_and_compute.hpp"


///
/// @param lastLoopEnd last file's loopEnd
/// @param newFlushNum
void mc_computation::parseSummary(size_t &lastLoopEnd, size_t &newFlushNum) {

//search for last_summary_observable.txt

    bool smrExists = false;
    std::regex smrRegex("last_summary_" + observable + ".txt");
    std::smatch matchSmr;
    for (const auto &entry: fs::directory_iterator(this->summary_folder)) {
        if (std::regex_search(entry.path().string(), matchSmr, smrRegex)) {
            smrExists = true;
            break;
        }


    }

    //search for last data file

    if (smrExists== false){
        newFlushNum=defaultFlushNum;
    }


}


///
/// @param pklDirectory directory containing pkl files
/// @param sortedFiles sorted pkl files by loop number
/// @param empty whether the data fils is empty
void mc_computation::sortOneDir(const std::string &pklDirectory,std::vector<std::string>& sortedFiles, bool &empty) {

    //search for pkl files
    std::regex pklRegex(".pkl");
    std::smatch matchPkl;
    std::vector<std::string> pklFilesAll;
    for (const auto &entry: fs::directory_iterator(pklDirectory)) {
        std::string entryStr = entry.path().string();
        if (std::regex_search(entryStr, matchPkl, pklRegex)) {
            pklFilesAll.push_back(entryStr);
        }
    }
    if (pklFilesAll.size()==0){

        empty=true;
        return;
    }

    std::vector<size_t> loopEndValsAll;
    std::regex lpEndRegex("loopEnd(\\d+)");
    std::smatch matchLpEnd;

    for (const std::string &name: pklFilesAll) {
        if (std::regex_search(name, matchLpEnd, lpEndRegex)) {
            loopEndValsAll.push_back(std::stoull(matchLpEnd.str(1)));
        }

    }

    std::vector<size_t> inds = argsort<size_t>(loopEndValsAll);


    for (const auto &i: inds) {
        sortedFiles.push_back(pklFilesAll[i]);
    }
    empty= false;
}



///
/// @param pklDirectory directory containing pkl files
/// @param lastPkl last pkl file
void mc_computation::searchLastDataInFolder(const std::string &pklDirectory, std::string & lastPkl) {
    std::vector<std::string> sortedFiles;
    bool empty;
    this->sortOneDir(pklDirectory, sortedFiles, empty);

    if (empty == true) {
        lastPkl = "";
        return;

    }

    size_t length = sortedFiles.size();
    lastPkl = sortedFiles[length - 1];


}


///
/// @param L_lastPklFile
/// @param y0_lastPklFile
/// @param z0_lastPklFile
/// @param y1_lastPklFile
void mc_computation::search_pkl(std::string& L_lastPklFile,std::string &y0_lastPklFile,
            std::string& z0_lastPklFile,std::string & y1_lastPklFile) {


    std::string L_pklStr;
    std::string y0_pklStr;
    std::string z0_pklStr;
    std::string y1_pklStr;

//search last pkl files
    this->searchLastDataInFolder(this->L_dataFolder, L_pklStr);
    this->searchLastDataInFolder(this->y0_dataFolder, y0_pklStr);
    this->searchLastDataInFolder(this->z0_dataFolder, z0_pklStr);
    this->searchLastDataInFolder(this->y1_dataFolder, y1_pklStr);


    //if there are no pkl files under each folder
    if(L_pklStr=="" and y0_pklStr=="" and z0_pklStr=="" and y1_pklStr==""){


        L_lastPklFile=L_pklStr;
        y0_lastPklFile=y0_pklStr;
        z0_lastPklFile=z0_pklStr;
        y1_lastPklFile=y1_pklStr;

        return;

    }


    //match loopLast

    size_t L_lpLast,y0_lpLast,z0_lpLast, y1_lpLast;
    std::regex lpEndRegex("loopEnd(\\d+)");
    std::smatch matchLpEnd;


    if(std::regex_search(L_pklStr,matchLpEnd,lpEndRegex)){
        L_lpLast=std::stoull(matchLpEnd.str(1));
    }

    if(std::regex_search(y0_pklStr,matchLpEnd,lpEndRegex)){
        y0_lpLast=std::stoull(matchLpEnd.str(1));
    }

    if(std::regex_search(z0_pklStr,matchLpEnd,lpEndRegex)){
        z0_lpLast=std::stoull(matchLpEnd.str(1));
    }

    if(std::regex_search(y1_pklStr,matchLpEnd,lpEndRegex)){
        y1_lpLast=std::stoull(matchLpEnd.str(1));
    }


    if(L_lpLast==y0_lpLast and y0_lpLast==z0_lpLast and z0_lpLast==y1_lpLast){

        L_lastPklFile=L_pklStr;
        y0_lastPklFile=y0_pklStr;
        z0_lastPklFile=z0_pklStr;
        y1_lastPklFile=y1_pklStr;

        return;
    }


    std::cerr<<"reading data error."<<std::endl;
    std::exit(2);


}


///
/// @param LCurr
/// @param y0Curr
/// @param z0Curr
/// @param y1Curr
void mc_computation::initialization(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr){
    double a = 10;
    LCurr=2*a;

    y0Curr=0.1*LCurr;

    z0Curr=0.2*LCurr;

    y1Curr=0.5*LCurr;


}

///
/// @param L
/// @param y0
/// @param z0
/// @param y1
/// @return
double mc_computation::f(const double &L,const double& y0, const double &z0, const double&y1){
    return this->beta * ((*potFuncPtr)(L,y0,z0,y1));


}