/// \file

#include <iostream>
#include "red_black.h"
#include "sequencial.h"

int main() {
    sequencial::poisson_gs(100);
    red_black::poisson_gs(100);
    return EXIT_SUCCESS;
}
