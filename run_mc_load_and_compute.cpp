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

    std::string x = "qwmn12";
    std::string y = "45iuqw";

    std::regex qwR("(\\d+)");
    std::smatch mc;

    if (std::regex_search(x, mc, qwR)) {
        std::cout << mc.str(1) << std::endl;
    }

    if (std::regex_search(y, mc, qwR)) {
        std::cout << mc.str(1) << std::endl;
    }

}