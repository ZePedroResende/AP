/// \file

#include <iostream>
#include "red_black.h"
#include "sequencial.h"

int main() {
    //    sequencial::poisson_gs(5);
    std::cout << std::endl;
    red_black::poisson_gs(5);
    return EXIT_SUCCESS;
}
