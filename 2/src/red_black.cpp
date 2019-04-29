#include "red_black.h"

void black(int max_i, int max_j, Matrice<double>& u, Matrice<double>& w) {
    for (int i = 1; i < max_i; i++)
        for (int j = 1; j < max_j; j++)
            if ((i + j) % 2 == 1)
                w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
}

void red(int max_i, int max_j, Matrice<double>& u, Matrice<double>& w) {
    for (int i = 1; i < max_i; i++)
        for (int j = 1; j < max_j; j++)
            if ((i + j) % 2 == 0)
                w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
}

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

        std::cout << "-----" << iter << "-----\n"
                  << w << "\n"
                  << u << "\n-----" << iter << "-----\n"
                  << std::endl;
        u = w;

        iter++;
    }

    std::cout << u << std::endl;
    std::cout << iter << std::endl;
}
