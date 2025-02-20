#include "red_black.h"

namespace {
void black(int max_i, int max_j, Matrice<double>& u, Matrice<double>& w) {
    int jstart;
    for (int i = 1; i < max_i; i++) {
        if (i % 2 == 1)
            jstart = 2;  // odd row
        else
            jstart = 1;  // even row
        for (int j = jstart; j < max_j; j += 2) {
            w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
        }
    }
}

void red(int max_i, int max_j, Matrice<double>& u, Matrice<double>& w) {
    int jstart;
    for (int i = 1; i < max_i; i++) {
        if (i % 2 == 1)
            jstart = 1;  // odd row
        else
            jstart = 2;  // even row
        for (int j = jstart; j < max_j; j += 2) {
            w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
        }
    }
}
}  // namespace

void red_black::poisson_gs(const int n) {
    double tol = pow((1.0 / n), 2);
    double diff = tol + 1;
    int iter = 0;

    Matrice<double> u(n, n);
    Matrice<double> w(n, n);
    Matrice<double> res(n, n);

    for (int i = 0; i < n - 1; i++) {
        u(0, i) = 100;
        w(0, i) = 100;
        u(i, 0) = 100;
        w(i, 0) = 100;
        u(i, n - 1) = 100;
        w(i, n - 1) = 100;
        u(n - 1, i) = 0;
        w(n - 1, i) = 0;
    }

    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < n - 1; j++) {
            u(i, j) = 50;
        }
    }

    while (diff > tol) {
        red(n - 1, n - 1, u, w);
        black(n - 1, n - 1, u, w);

        diff = fabs((w - u).max());
        u = w;
        iter++;
    }

    //    std::cout << u << std::endl;
    //   std::cout << iter << std::endl;
}
