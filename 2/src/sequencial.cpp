#include "sequencial.h"

void sequencial::poisson_gs(const int n) {
    int tol = (1 / n) ^ 2;
    int diff = tol + 1;
    int iter = 0;

    Matrice<int> u(n, n);
    Matrice<int> w(n, n);

    for (int i = 0; i < n; i++) {
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
        for (int i = 1; i < n - 1; i++)
            for (int j = 1; j < n - 1; j++)
                w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4;

        diff = (w - u).max();
        u = w;
        iter++;
    }

    std::cout << u << std::endl;
    //    std::cout << iter << std::endl;
}
