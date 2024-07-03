#include "mc_load_and_compute.hpp"



void mc_computation::parseSummary(size_t &newFlushNum,bool&moreDataPointsAfterEq) {

//search for last_summary_observable.txt

    bool smrExists = false;
    std::string fileName;
    std::regex smrRegex("summaryFile_" + observable + ".txt");
    std::smatch matchSmr;
    size_t lag = 0;
    size_t newDataPointNum = 0;
    for (const auto &entry: fs::directory_iterator(this->summary_folder)) {
        if (std::regex_search(entry.path().string(), matchSmr, smrRegex)) {
            smrExists = true;
            fileName = entry.path().string();
            break;
        }


    }

    //search for summary file

    if (smrExists == false) {
        newFlushNum = defaultFlushNum;
        moreDataPointsAfterEq = false;
        return;
    }

    //parse summary file

    std::regex eqRegex("equilibrium");
    std::regex errRegex("error");
    std::regex lagRegex("lag=(\\d+)");
    std::regex newDataNum("newDataPointNum=(\\d+)");
    std::regex continueRegex("continue");
    std::regex highCorrRegex("high");

    std::smatch matchEq;
    std::smatch matchErr;
    std::smatch matchLag;
    std::smatch matchNewNum;
    std::smatch matchContinue;
    std::smatch matchHigh;


    std::ifstream ifs(fileName);
    if (!ifs.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        std::exit(1);
    }

    std::string line;
    while (std::getline(ifs, line)) {
        if (std::regex_search(line, matchErr, errRegex)) {
            std::cerr << "error in previous computation, please re-run." << std::endl;
            std::exit(10);
        }//end of matching error

        if (std::regex_search(line, matchEq, eqRegex)) {
            continue;

        }// end of matching equilibrium

        if (std::regex_search(line, matchContinue, continueRegex)) {
            newFlushNum = this->defaultFlushNum;
            moreDataPointsAfterEq = false;
            return;
        }//end of matching continue

        if (std::regex_search(line, matchHigh, highCorrRegex)) {
            newFlushNum = this->defaultFlushNum;
            moreDataPointsAfterEq = false;
            return;
        }//end of matching high

        if (std::regex_search(line, matchLag, lagRegex)) {
            lag = std::stoull(matchLag.str(1));
        }//end of  matching lag

        if (std::regex_search(line, matchNewNum, newDataNum)) {
            newDataPointNum = std::stoull(matchNewNum.str(1));
        }//end of matching newDataPointNum

    }


    if (newDataPointNum == 0) {
        newFlushNum = 0;
        moreDataPointsAfterEq = true;
        return;
    } else {
        size_t totalNum = lag * newDataPointNum;
        newFlushNum = static_cast<size_t>(std::ceil(
                static_cast<double>(totalNum) / (static_cast<double>(loopToWrite))));
        moreDataPointsAfterEq = true;
        return;
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
void mc_computation::initializeParams(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr){
    double a = 10;
    LCurr=2*a;

    y0Curr=0.1*LCurr;

    z0Curr=0.2*LCurr;

    y1Curr=0.5*LCurr;


}




///
/// @param LInit initial value of L
/// @param y0Init initial value of y0
/// @param z0Init initial value of z0
/// @param y1Init initial value of y1
void mc_computation::execute_mc(const double& LInit,const double &y0Init, const double &z0Init, const double& y1Init, const size_t & loopInit, const size_t & flushNum) {
    double LCurr=LInit;
    double y0Curr=y0Init;
    double z0Curr=z0Init;
    double y1Curr=y1Init;
    double UCurr=(*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);

    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    size_t loopStart = loopInit;
    for(size_t fls=0;fls<flushNum;fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};

        for (size_t j = 0; j < loopToWrite; j++) {
            //propose a move
            double LNext;
            double y0Next;
            double z0Next;
            double y1Next;
            this->proposal(LCurr, y0Curr, z0Curr, y1Curr, LNext, y0Next, z0Next, y1Next);
            double UNext;
            double r = acceptanceRatio(LCurr, y0Curr, z0Curr, y1Curr, UCurr, LNext, y0Next, z0Next, y1Next, UNext);
            double u = distUnif01(e2);
            if (u <= r) {
                LCurr = LNext;
                y0Curr = y0Next;
                z0Curr = z0Next;
                y1Curr = y1Next;
                UCurr = UNext;

            }//end of accept-reject

            U_ptr[j] = UCurr;
            L_ptr[j] = LCurr;
            y0_ptr[j] = y0Curr;
            z0_ptr[j] = z0Curr;
            y1_ptr[j] = y1Curr;
        }//end for loop


        size_t loopEnd = loopStart + (fls+1)*loopToWrite - 1;

        std::string fileNameMiddle = "loopStart" + std::to_string(loopStart) + "loopEnd" + std::to_string(loopEnd);
        //save U_ptr
        std::string out_UPickleFileName = this->U_dataFolder + "/" + fileNameMiddle + ".U.pkl";
        save_array_to_pickle(U_ptr, loopToWrite, out_UPickleFileName);

        //save L_ptr
        std::string out_LPickleFileName = L_dataFolder + "/" + fileNameMiddle + ".L.pkl";
        save_array_to_pickle(L_ptr, loopToWrite, out_LPickleFileName);

        //save y0_ptr
        std::string out_y0PickleFileName = y0_dataFolder + "/" + fileNameMiddle + ".y0.pkl";
        save_array_to_pickle(y0_ptr, loopToWrite, out_y0PickleFileName);

        //save z0_ptr
        std::string out_z0PickleFileName = z0_dataFolder + "/" + fileNameMiddle + ".z0.pkl";
        save_array_to_pickle(z0_ptr, loopToWrite, out_z0PickleFileName);

        //save y1_ptr
        std::string out_y1PickleFileName = y1_dataFolder + "/" + fileNameMiddle + ".y1.pkl";
        save_array_to_pickle(y1_ptr, loopToWrite, out_y1PickleFileName);

        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout<<"loop "+std::to_string(loopStart)+" to loop "+std::to_string(loopEnd)+": "<<elapsed_secondsAll.count()/3600.0 << " h" << std::endl;

        loopStart=loopEnd+1;

    }//end flush for loop

}


///
/// @param LCurr current value of L
/// @param y0Curr current value of y0
/// @param z0Curr current value of z0
/// @param y1Curr current value of y1
/// @param LNext  next value of L
/// @param y0Next next value of y0
/// @param z0Next next value of z0
/// @param y1Next next value of y1
void mc_computation::proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
              double & LNext, double & y0Next, double & z0Next, double & y1Next){

//next L
LNext= generate_nearby_normal(LCurr,this->h);

//next y0
y0Next= generate_nearby_normal(y0Curr,this->h);

//next z0
z0Next= generate_nearby_normal(z0Curr,this->h);

//next y1
y1Next= generate_nearby_normal(y1Curr,this->h);

}

double mc_computation::acceptanceRatio(const double &LCurr,const double &y0Curr,
                                              const double &z0Curr, const double& y1Curr,const double& UCurr,
                                              const double &LNext, const double& y0Next,
                                              const double & z0Next, const double & y1Next,
                                              double &UNext){

    UNext=((*potFuncPtr)(LNext,y0Next,z0Next,y1Next));

    double numerator = -this->beta*UNext;

    double denominator=-this->beta*UCurr;

    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}

void mc_computation::loadParams(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr,
                const std::string & L_lastPklFile, const std::string &y0_lastPklFile,
                const std::string &z0_lastPklFile, const std::string & y1_lastPklFile){

    //load L
    size_t last_L_size;
    std::shared_ptr<double[]> last_L_ptr=load_array_from_pickle(last_L_size,L_lastPklFile);
    LCurr=last_L_ptr[last_L_size-1];

    //load y0
    size_t  last_y0_size;
    std::shared_ptr<double[]> last_y0_ptr=load_array_from_pickle(last_y0_size,y0_lastPklFile);
    y0Curr=last_y0_ptr[last_y0_size-1];

    //load z0
    size_t last_z0_size;
    std::shared_ptr<double[]> last_z0_ptr= load_array_from_pickle(last_z0_size,z0_lastPklFile);
    z0Curr=last_z0_ptr[last_z0_size-1];

    //load y1
    size_t  last_y1_size;
    std::shared_ptr<double[]> last_y1_ptr= load_array_from_pickle(last_y1_size,y1_lastPklFile);
    y1Curr=last_y1_ptr[last_y1_size-1];




}

void mc_computation::strategy(double &LInit, double &y0Init, double &z0Init, double& y1Init, size_t& flushNum,size_t & loopInit,bool&moreDataPointsAfterEq) {

    std::string L_lastPklFile;
    std::string y0_lastPklFile;
    std::string z0_lastPklFile;
    std::string y1_lastPklFile;

    search_pkl(L_lastPklFile, y0_lastPklFile, z0_lastPklFile, y1_lastPklFile);
//    std::cout<<"L_lastPklFile="<<L_lastPklFile<<std::endl;
    if (L_lastPklFile == "" and y0_lastPklFile == "" and z0_lastPklFile == ""
        and y1_lastPklFile == "") {
        this->initializeParams(LInit, y0Init, z0Init, y1Init);
        flushNum = defaultFlushNum;
        loopInit = 0;
        moreDataPointsAfterEq= false;
        return;
    } else {
        this->parseSummary(flushNum,moreDataPointsAfterEq);

        this->loadParams(LInit, y0Init, z0Init, y1Init, L_lastPklFile,
                         y0_lastPklFile, z0_lastPklFile, y1_lastPklFile);


        std::regex lpEndRegex("loopEnd(\\d+)");
        std::smatch matchLpEnd;
        if (std::regex_search(L_lastPklFile, matchLpEnd, lpEndRegex)) {

            loopInit = std::stoull(matchLpEnd.str(1)) + 1;
        }


        return;

    }


}


void mc_computation::load_init_run() {
    double LInit;
    double y0Init;
    double z0Init;
    double y1Init;
    size_t flushNum;
    size_t loopInit;
    bool moreDataPointsAfterEq;
    this->strategy(LInit, y0Init, z0Init, y1Init, flushNum, loopInit, moreDataPointsAfterEq);


    std::cout << "LInit=" << LInit << std::endl;
    std::cout << "y0Init=" << y0Init << std::endl;
    std::cout << "z0Init=" << z0Init << std::endl;
    std::cout << "y1Init=" << y1Init << std::endl;
    std::cout << "flushNum=" << flushNum << std::endl;
    std::cout << "loopInit=" << loopInit << std::endl;

    execute_mc(LInit, y0Init, z0Init, y1Init, loopInit, flushNum);

    std::string smrFile = summary_folder + "/summaryFile_" + observable + ".txt";
    if (moreDataPointsAfterEq == true) {
        std::cout << "mc executed for " << flushNum << " more flushes." << std::endl;
        std::ofstream outFile(smrFile, std::ios::app);
        if (outFile.is_open()) {
            // Write the text to the file
            outFile << "mc executed for " << flushNum << " more flushes." << std::endl;

            // Close the file
            outFile.close();
        } else {
            std::cerr << "Unable to open file." << std::endl;
        }

    }//end of writing to summary file


}