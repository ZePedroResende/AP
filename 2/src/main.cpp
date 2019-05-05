/// \file

#include <chrono>
#include <iostream>
#include "block_independent.h"
#include "independent.h"
#include "red_black.h"
#include "sequencial.h"

int main(int argc, char* argv[]) {
    const int n = 100;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    sequencial::poisson_gs(n);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
              << std::endl;
    begin = std::chrono::steady_clock::now();
    red_black::poisson_gs(n);
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
              << std::endl;
    begin = std::chrono::steady_clock::now();
    independent::poisson_gs(n);
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
              << std::endl;
    begin = std::chrono::steady_clock::now();
    block_independent::poisson_gs(argc, argv, n);
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
              << std::endl;

    return EXIT_SUCCESS;
}
