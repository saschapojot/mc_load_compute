//
// Created by polya on 6/29/24.
//

#ifndef MC_LOAD_COMPUTE_MC_LOAD_AND_COMPUTE_HPP
#define MC_LOAD_COMPUTE_MC_LOAD_AND_COMPUTE_HPP
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/object/pickle_support.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <chrono>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <math.h>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
const auto PI=M_PI;



class potentialFunction {
//base class for potential function
public:
    virtual double operator()(const double&L,const double &y0, const double &z0, const double& y1) const = 0;
    virtual std::string getParamStr() const = 0; // Pure virtual function
    virtual ~ potentialFunction() {};
};


class quadratic : public potentialFunction {
public:
    quadratic(const std::string &row) : potentialFunction() {
//        this->a1 = 1;
//        this->a2 = 1.5;
//        this-> c1 = 50;
//        this->c2 = 80;
        this->rowName = row;
        this->parseParams();


        this->paramStr = rowName;
        std::cout << "a1=" << a1 << ", a2=" << a2 << ", c1=" << c1 << ", c2=" << c2 << std::endl;


    }

public:
    std::string getParamStr() const override {
        return this->paramStr;
    }

public:
    double operator()(const double &L, const double &y0, const double &z0, const double &y1) const override {
        double val = c1 * std::pow(y0 - a1, 2) + c2 * std::pow(z0 - a2, 2)
                     + c1 * std::pow(y1 - a1, 2) + c2 * std::pow(-y0 - z0 - y1 + L - a2, 2);
//        std::cout<<"val="<<val<<std::endl;
        return val;

    }


    void parseParams() {

        std::regex rowRegex("(row\\d++)");
        std::regex a1Regex("a1\\s*=\\s*([+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?)");
        std::regex a2Regex("a2\\s*=\\s*([+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?)");
        std::regex c1Regex("c1\\s*=\\s*([+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?)");
        std::regex c2Regex("c2\\s*=\\s*([+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?)");

        std::smatch matchRow;
        std::smatch match_a1;
        std::smatch match_a2;
        std::smatch match_c1;
        std::smatch match_c2;

        std::ifstream file(inParamsFile);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << inParamsFile << std::endl;
            std::exit(1);
        }
        bool rowFound = false;
        std::string line;
        while (std::getline(file, line)) {
            // Process the line (e.g., print it)
//            std::cout << line << std::endl;
            //match row
            if (std::regex_search(line, matchRow, rowRegex)) {
                if (rowName == matchRow.str(1)) {
                    rowFound = true;
                }

                //match params
                //match a1
                if (std::regex_search(line, match_a1, a1Regex)) {
                    this->a1 = std::stod(match_a1.str(1));
                } else {
                    std::cerr << "a1 missing." << std::endl;
                    std::exit(2);
                }

                //match a2
                if (std::regex_search(line, match_a2, a2Regex)) {
                    this->a2 = std::stod(match_a2.str(1));
                } else {
                    std::cerr << "a2 missing." << std::endl;
                    std::exit(2);
                }
                //match c1
                if (std::regex_search(line, match_c1, c1Regex)) {
                    this->c1 = std::stod(match_c1.str(1));
                } else {
                    std::cerr << "c1 missing." << std::endl;
                    std::exit(2);
                }

                //match c2
                if (std::regex_search(line, match_c2, c2Regex)) {
                    this->c2 = std::stod(match_c2.str(1));
                } else {
                    std::cerr << "c2 missing." << std::endl;
                    std::exit(2);
                }
            }
            break;

        }//end of while

        if (rowFound == false) {
            std::cerr << rowName + " not found." << std::endl;
            std::exit(3);
        }
    }// end of parseParams()




public:

    double a1;
    double a2;
    double c1;
    double c2;
    std::string paramStr;
    std::string rowName;
    std::string inParamsFile = "./inputData/quadratic/quadraticCoeffs.txt";

};


class mc_computation{
public:
    mc_computation(const double& temperature, const std::shared_ptr<potentialFunction> &funcPtr,const  std::string &observableName){

        this->T = temperature;
        this->beta = 1 / T;
        this->potFuncPtr = funcPtr;

        double stepForT1 = 0.1;
        this->h = stepForT1 * T>0.2? 0.2:stepForT1 * T;//stepSize;
        std::cout << "h=" << h << std::endl;




    ///create directories if not exist

     this->funcName = this->demangle(typeid(*potFuncPtr).name());
     this->observable=observableName;


     this->TFolder=data_root+"/"+funcName+"/"+potFuncPtr->getParamStr()+"/T"+std::to_string(T)+"/";

        //create directory TFolder if it does not exist
        if (!fs::is_directory(TFolder) || !fs::exists(TFolder)) {
            fs::create_directories(TFolder);
        }


        this->data_folder=TFolder+"/data_files/";
        //create directory data_folder if it does not exist
        if (!fs::is_directory(data_folder) || !fs::exists(data_folder)) {
            fs::create_directories(data_folder);
        }

        this->U_dataFolder=data_folder+"/U_AllPickle/";
        //create directory U_dataFolder if it does not exist
        if (!fs::is_directory(U_dataFolder) || !fs::exists(U_dataFolder)) {
            fs::create_directories(U_dataFolder);
        }


        this->L_dataFolder=data_folder+"/L_AllPickle/";
        //create directory L_dataFolder if it does not exist
        if (!fs::is_directory(L_dataFolder) || !fs::exists(L_dataFolder)) {
            fs::create_directories(L_dataFolder);
        }


        this->y0_dataFolder=data_folder+"/y0_AllPickle/";
        //create directory y0_dataFolder if it does not exist
        if (!fs::is_directory(y0_dataFolder) || !fs::exists(y0_dataFolder)) {
            fs::create_directories(y0_dataFolder);
        }

        this->z0_dataFolder=data_folder+"/z0_AllPickle/";
        //create directory z0_dataFolder if it does not exist
        if (!fs::is_directory(z0_dataFolder) || !fs::exists(z0_dataFolder)) {
            fs::create_directories(z0_dataFolder);
        }


        this->y1_dataFolder=data_folder+"/y1_AllPickle/";
        //create directory y1_dataFolder if it does not exist
        if (!fs::is_directory(y1_dataFolder) || !fs::exists(y1_dataFolder)) {
            fs::create_directories(y1_dataFolder);
        }



        this->summary_folder=TFolder+"/summary_"+observable+"/";
        //create directory summary_folder if it does not exist
        if (!fs::is_directory(summary_folder) || !fs::exists(summary_folder)) {
            fs::create_directories(summary_folder);
        }

//        this->summary_history_folder=summary_folder+"/summary_history/";
//        //create directory summary_history_folder if it does not exist
//        if (!fs::is_directory(summary_history_folder) || !fs::exists(summary_history_folder)) {
//            fs::create_directories(summary_history_folder);
//        }


        ///allocate arrays
        try {
            U_ptr = std::shared_ptr<double[]>(new double[loopToWrite],
                                              std::default_delete<double[]>());
            L_ptr = std::shared_ptr<double[]>(new double[loopToWrite],
                                              std::default_delete<double[]>());
            y0_ptr = std::shared_ptr<double[]>(new double[loopToWrite],
                                               std::default_delete<double[]>());
            z0_ptr = std::shared_ptr<double[]>(new double[loopToWrite],
                                               std::default_delete<double[]>());
            y1_ptr = std::shared_ptr<double[]>(new double[loopToWrite],
                                               std::default_delete<double[]>());


        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }

    }



public:

    //generate the name of the potential function
    std::string demangle(const char *name) {
        int status = -1;
        char *demangled = abi::__cxa_demangle(name, NULL, NULL, &status);
        std::string result(name);
        if (status == 0) {
            result = demangled;
        }
        std::free(demangled);
        return result;
    }

    void load_init_run();

    ///
    /// @param LInit initial value of L
    /// @param y0Init initial value of y0
    /// @param z0Init initial value of z0
    /// @param y1Init initial value of y1
    /// @param flushNum files to write
    void execute_mc(const double& LInit,const double &y0Init, const double &z0Init, const double& y1Init, const size_t & loopInit, const size_t & flushNum);


    ///
    /// @param LCurr current value of L
    /// @param y0Curr current value of y0
    /// @param z0Curr current value of z0
    /// @param y1Curr current value of y1
    /// @param LNext  next value of L
    /// @param y0Next next value of y0
    /// @param z0Next next value of z0
    /// @param y1Next next value of y1
    void proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
                  double & LNext, double & y0Next, double & z0Next, double & y1Next);


    ///
    /// @param x
    /// @param sigma
    /// @return a value around x, from a  normal distribution
    static double generate_nearby_normal(const double & x, const double &sigma){
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(x, sigma); // Normal distribution with mean rCurr and standard deviation sigma

        double xNext = d(gen);


        return xNext;


    }

    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    /// @param LNext
    /// @param y0Next
    /// @param z0Next
    /// @param y1Next
    /// @param UNext
    /// @return
    double acceptanceRatio(const double &LCurr,const double &y0Curr,
                           const double &z0Curr, const double& y1Curr,const double& UCurr,
                           const double &LNext, const double& y0Next,
                           const double & z0Next, const double & y1Next,
                           double &UNext);

    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    void initializeParams(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr);

    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    /// @param L_lastPklFile
    /// @param y0_lastPklFile
    /// @param z0_lastPklFile
    /// @param y1_lastPklFile
    void loadParams(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr,
                    const std::string & L_lastPklFile, const std::string &y0_lastPklFile,
                    const std::string &z0_lastPklFile, const std::string & y1_lastPklFile);

    ///
    /// @param L_lastPklFile
    /// @param y0_lastPklFile
    /// @param z0_lastPklFile
    /// @param y1_lastPklFile
    void search_pkl(std::string& L_lastPklFile,std::string &y0_lastPklFile,
                std::string& z0_lastPklFile,std::string & y1_lastPklFile);


    void strategy(double &LInit, double &y0Init, double &z0Init, double& y1Init, size_t& flushNum,size_t & loopInit);

    ///
    /// @param lastLoopEnd last file's loopEnd
    /// @param newFlushNum
    void parseSummary(size_t &lastLoopEnd, size_t &newFlushNum);


    ///
    /// @param pklDirectory directory containing pkl files
    /// @param lastPkl last pkl file
    void searchLastDataInFolder(const std::string &pklDirectory, std::string & lastPkl);


    ///
    /// @param pklDirectory directory containing pkl files
    /// @param sortedFiles sorted pkl files by loop number
    /// @param empty whether the data fils is empty
   void sortOneDir(const std::string &pklDirectory,
                   std::vector<std::string>& sortedFiles, bool &empty);



    template<class T>
    static std::vector<size_t> argsort(const std::vector<T> &v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] <= v[i2]; });
        return idx;
    }

    //read and write to pickle
    static void save_array_to_pickle(std::shared_ptr<double[]> &ptr, std::size_t size, const std::string& filename) {
        using namespace boost::python;
        try {
            Py_Initialize();  // Initialize the Python interpreter
            if (!Py_IsInitialized()) {
                throw std::runtime_error("Failed to initialize Python interpreter");
            }

            // Debug output
            std::cout << "Python interpreter initialized successfully." << std::endl;

            // Import the pickle module
            object pickle = import("pickle");
            object pickle_dumps = pickle.attr("dumps");

            // Create a Python list from the C++ array
            list py_list;
            for (std::size_t i = 0; i < size; ++i) {
                py_list.append(ptr[i]);
            }

            // Serialize the list using pickle.dumps
            object serialized_array = pickle_dumps(py_list, 2);  // Use protocol 2 for binary compatibility

            // Extract the serialized data as a string
            std::string serialized_str = extract<std::string>(serialized_array);

            // Write the serialized data to a file
            std::ofstream file(filename, std::ios::binary);
            if (!file) {
                throw std::runtime_error("Failed to open file for writing: " + filename);
            }
            file.write(serialized_str.data(), serialized_str.size());
            file.close();

            // Debug output
            std::cout << "Array serialized and written to file successfully." << std::endl;
        } catch (const error_already_set&) {
            PyErr_Print();
            std::cerr << "Boost.Python error occurred while saving array to pickle file." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }

        if (Py_IsInitialized()) {
            Py_Finalize();  // Finalize the Python interpreter
        }
    }


    static std::shared_ptr<double[]> load_array_from_pickle(std::size_t& size, const std::string& filename) {
        using namespace boost::python;
        std::shared_ptr<double[]> ptr;

        try {
            Py_Initialize();  // Initialize the Python interpreter
            if (!Py_IsInitialized()) {
                throw std::runtime_error("Failed to initialize Python interpreter");
            }

            // Debug output
            std::cout << "Python interpreter initialized successfully." << std::endl;

            // Import the pickle module
            object pickle = import("pickle");
            object pickle_loads = pickle.attr("loads");

            // Read the serialized data from the file
            std::ifstream file(filename, std::ios::binary);
            if (!file) {
                throw std::runtime_error("Failed to open file for reading: " + filename);
            }
            std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
            file.close();

            // Create a Python bytes object from the buffer
            object py_bytes = object(handle<>(PyBytes_FromStringAndSize(buffer.data(), buffer.size())));

            // Deserialize the data using pickle.loads
            object py_list = pickle_loads(py_bytes);

            // Get the size of the list
            size = len(py_list);
            ptr = std::shared_ptr<double[]>(new double[size]);

            // Convert the Python list to a C++ array
            for (std::size_t i = 0; i < size; ++i) {
                ptr[i] = extract<double>(py_list[i]);
            }

            // Debug output
            std::cout << "Array deserialized and loaded from file successfully." << std::endl;
        } catch (const error_already_set&) {
            PyErr_Print();
            std::cerr << "Boost.Python error occurred while loading array from pickle file." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }

        if (Py_IsInitialized()) {
            Py_Finalize();  // Finalize the Python interpreter
        }

        return ptr;
    }
public:
    double T;// temperature
    double beta;

    static const size_t loopToWrite=1000000;
    double h;// step size
    std::shared_ptr<potentialFunction> potFuncPtr;

    size_t defaultFlushNum=5;

    std::string data_folder;
    std::string U_dataFolder;
    std::string L_dataFolder;
    std::string y0_dataFolder;
    std::string z0_dataFolder;
    std::string y1_dataFolder;

    std::string  summary_folder;

    std::string data_root="./dataAll/";
    std::string funcName;

    std::string observable;
    std::string TFolder;
    std::string summary_history_folder;

    std::shared_ptr<double[]> U_ptr;
    std::shared_ptr<double[]> L_ptr;
    std::shared_ptr<double[]> y0_ptr;
    std::shared_ptr<double[]> z0_ptr;
    std::shared_ptr<double[]> y1_ptr;




};


#endif //MC_LOAD_COMPUTE_MC_LOAD_AND_COMPUTE_HPP
