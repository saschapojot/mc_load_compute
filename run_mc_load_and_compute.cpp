#include "mc_subroutine/mc_load_and_compute.hpp"


int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
    double T;
    try {
        T = std::stod(argv[1]); // Attempt to convert string to double

    } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
    } catch (const std::out_of_range &e) {
        std::cerr << "Out of range: " << e.what() << std::endl;
    }

    if (T <= 0) {
        std::cerr << "T must be >0" << std::endl;
        std::exit(1);
    }

    std::string observableName = std::string(argv[2]);

    auto mcObj = mc_computation(T, std::make_shared<quadratic>(), observableName);
    mcObj.load_init_run();
}