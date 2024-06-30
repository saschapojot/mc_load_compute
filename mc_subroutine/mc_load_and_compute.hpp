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
    quadratic() : potentialFunction() {
        this->a1 = 1;
        this->a2 = 1.5;
        this-> c1 = 50;
        this->c2 = 80;

        this->paramStr = "a1_" + std::to_string(a1) + "a2_" + std::to_string(a2) + "c1_" + std::to_string(c1)
                   + "c2_" + std::to_string(c2);

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

public:

    double a1;
    double a2;
    double c1;
    double c2;
    std::string paramStr;

};


class mc_computation{
public:
    mc_computation(const double& temperature, const std::shared_ptr<potentialFunction> &funcPtr,const  std::string &observableName){

        this->T = temperature;
        this->beta = 1 / T;
        this->potFuncPtr = funcPtr;






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

        this->summary_history_folder=summary_folder+"/summary_history/";
        //create directory summary_history_folder if it does not exist
        if (!fs::is_directory(summary_history_folder) || !fs::exists(summary_history_folder)) {
            fs::create_directories(summary_history_folder);
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

    void execute_mc();


    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    void initialization(double &LCurr, double &y0Curr, double &z0Curr, double &y1Curr);



    ///
    /// @param L
    /// @param y0
    /// @param z0
    /// @param y1
    /// @return
    double f(const double &L,const double& y0, const double &z0, const double&y1);
    ///
    /// @param L_lastPklFile
    /// @param y0_lastPklFile
    /// @param z0_lastPklFile
    /// @param y1_lastPklFile
    void search_pkl(std::string& L_lastPklFile,std::string &y0_lastPklFile,
                std::string& z0_lastPklFile,std::string & y1_lastPklFile);

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
            object serialized_array = pickle_dumps(py_list);

            // Extract the serialized data as a string
            std::string serialized_str = extract<std::string>(serialized_array);

            // Write the serialized data to a file
            std::ofstream file(filename, std::ios::binary);
            if (!file) {
                throw std::runtime_error("Failed to open file for writing");
            }
            file.write(serialized_str.data(), serialized_str.size());
            file.close();

            // Debug output
            std::cout << "Array serialized and written to file successfully." << std::endl;
        } catch (const error_already_set&) {
            PyErr_Print();
            std::cerr << "Boost.Python error occurred." << std::endl;
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
                throw std::runtime_error("Failed to open file for reading");
            }
            std::string serialized_str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
            file.close();

            // Deserialize the data using pickle.loads
            object py_list = pickle_loads(serialized_str);

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
            std::cerr << "Boost.Python error occurred." << std::endl;
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

    static const size_t loopToWrite=100000;
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




};


#endif //MC_LOAD_COMPUTE_MC_LOAD_AND_COMPUTE_HPP
