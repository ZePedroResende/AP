#include "independent.h"

namespace {
double black(const int max_i, const int max_j, Matrice<double>& u, Matrice<double>& w) {
    int jstart;
    double max, aux;
    for (int i = 1; i < max_i; i++) {
        if (i % 2 == 1)
            jstart = 2;  // odd row
        else
            jstart = 1;  // even row
        for (int j = jstart; j < max_j; j += 2) {
            w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
            aux = fabs(w(i, j) - u(i, j));
            max = max > aux ? max : aux;
        }
    }
    return max;
}

double red(const int max_i, const int max_j, Matrice<double>& u, Matrice<double>& w) {
    int jstart;
    double max, aux;
    for (int i = 1; i < max_i; i++) {
        if (i % 2 == 1)
            jstart = 1;  // odd row
        else
            jstart = 2;  // even row
        for (int j = jstart; j < max_j; j += 2) {
            w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
            aux = fabs(w(i, j) - u(i, j));
            max = max > aux ? max : aux;
        }
    }
    return max;
}
}  // namespace

void independent::poisson_gs(const int n) {
    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    double tol = pow((1.0 / n), 2);
    double diff = tol + 1;
    double a, b;
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
        a = red(n - 1, n - 1, u, w);
        b = black(n - 1, n - 1, u, w);
        diff = a > b ? a : b;
        u = w;
        iter++;
    }
    /*
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
              << std::endl;
              */
    //    std::cout << u << std::endl;
    //   std::cout << iter << std::endl;
}
