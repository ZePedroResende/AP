#ifndef BLOCK_INDEPENDENT_H
#define BLOCK_INDEPENDENT_H
#include <mpi.h>
#include <chrono>
#include <cmath>
#include <cstdlib> /* abs */
#include <iostream>
#include "matrice.h"

namespace block_independent {
void poisson_gs(int argc, char* argv[], const int n);
}
#endif  // BLOCK_INDEPENDENT_H
