#include "block_independent.h"

namespace {
double black(const int I, const int I_MAX, const int J, const int J_MAX,
             int rank, Matrice<double>& u, Matrice<double>& w) {
  int jstart;
  double max = -1;
  double aux;
  //  std::cout << w << u << std::endl;

  for (int i = (I + 1); i < I_MAX; i++) {
    if (i % 2 == 1)
      jstart = 2 + J;  // odd row
    else
      jstart = 1 + J;  // even row
    for (int j = jstart; j < J_MAX; j += 2) {
      w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
      aux = fabs(w(i, j) - u(i, j));
      max = max > aux ? max : aux;
    }
  }
  return max;
}

double red(const int I, const int I_MAX, const int J, const int J_MAX,
           const int rank, Matrice<double>& u, Matrice<double>& w) {
  int jstart;
  double max = -1;
  double aux;
  for (int i = (I + 1); i < I_MAX; i++) {
    if (i % 2 == 1)
      jstart = 1 + J;  // odd row
    else
      jstart = 2 + J;  // even row
    for (int j = jstart; j < J_MAX; j += 2) {
      w(i, j) = (w(i - 1, j) + w(i, j - 1) + u(i, j + 1) + u(i + 1, j)) / 4.0;
      aux = fabs(w(i, j) - u(i, j));
      max = max > aux ? max : aux;
    }
  }
  return max;
}
}  // namespace

void init_matrices(Matrice<double>& u, Matrice<double>& w, int n) {
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
}

void comunicate_collumns(Matrice<double>& w, int I, int I_MAX, int J, int J_MAX,
                         int I_, int I_MAX_, int J_, int J_MAX_, int rank,
                         int size) {
  std::vector<double> v;
  MPI_Status status;

  v = w.get_vector(I, I_MAX, J, J_MAX);
  MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
  v.clear();
  v.resize(size);
  MPI_Recv(&v[0], size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
  w.update_with_vector(I_, I_MAX_, J_, J_MAX_, v);
}

void update_collumns(Matrice<double>& w, int I, int I_MAX, int J, int J_MAX,
                     int rank, int n, int nb) {
  int i, i_max, j, j_max;
  int i_, i_max_, j_, j_max_;
  i = I == 0 ? I : I + 1;
  i_max = I_MAX == n - 1 ? I_MAX : I_MAX - 1;
  j = J == 0 ? J : J + 1;
  j_max = J_MAX == n - 1 ? J_MAX : J_MAX - 1;

  i_ = I == 0 ? I : I + 1;
  i_max_ = I_MAX == n - 1 ? I_MAX : I_MAX - 1;
  j_ = J == 0 ? J : J + 1;
  j_max_ = J_MAX == n - 1 ? J_MAX : J_MAX - 1;

  if (I != 0) {
    comunicate_collumns(w, I + 1, I + 2, j, j_max, I, I + 1, j_, j_max_,
                        rank - nb, j_max_ - j_);
  }
  if (I_MAX != n - 1) {
    comunicate_collumns(w, I_MAX - 2, I_MAX - 1, j, j_max, I_MAX - 1, I_MAX, j_,
                        j_max_, rank + nb, j_max_ - j_);
  }
  if (J != 0) {
    comunicate_collumns(w, i, i_max, J + 1, J + 2, i_, i_max_, J, J + 1,
                        rank - 1, i_max_ - i_);
  }
  if (J_MAX != n - 1) {
    comunicate_collumns(w, i, i_max, J_MAX - 2, J_MAX - 1, i_, i_max_,
                        J_MAX - 1, J_MAX, rank + 1, i_max_ - i_);
  }
}

Matrice<double>& mpi_logic(Matrice<double>& u, Matrice<double>& w,
                           Matrice<double>& res, int n, int rank, int NP) {
  double tol = pow((1.0 / n), 2);
  double diff = tol + 1;
  // std::cout << "diff= " << diff << std::endl;
  double a, b;
  int iter = 0;
  // mpi
  int P = NP / 2;
  int X = std::floor(sqrt(P));
  int k = pow(X, 2);
  int msg;
  MPI_Status status;

  std::vector<double> v;
  int nb = std::floor(sqrt(k));
  int B = std::floor(n / sqrt(k));
  int J = B * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
              ? 0
              : B * (rank - (nb * (int)std::floor(rank / nb))) - 1;
  int J_MAX = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                  ? B * (rank - (nb * (rank / nb)) + 1)
                  : n - 1;
  int rank_j = rank < k ? rank : rank - k;
  int I = std::floor(rank_j / nb) * B - 1 == -1 ? 0 : (rank_j / nb) * B - 1;
  int I_MAX = ((rank_j / nb) + 1) < nb ? ((rank_j / nb) + 1) * B : n - 1;

  if (rank == 0) {
    for (;;) {
      a = red(I, I_MAX, J, J_MAX, rank, u, w);
      MPI_Recv(&b, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
      a = b < a ? a : b;

      v.clear();
      v.resize((I_MAX - I) * (J_MAX - J));
      MPI_Recv(&v[0], v.size(), MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
      w.update_black_vector(I, I_MAX, J, J_MAX, v);

      for (int l = 1; l < k; l++) {
        MPI_Recv(&b, 1, MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);
        a = a < b ? b : a;
      }

      diff = a;
      msg = diff <= tol;
      MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (msg) {
        for (int l = 1; l < k; l++) {
          int J_ = B * (l - (nb * (int)std::floor(l / nb))) - 1 == -1
                       ? 0
                       : B * (l - (nb * (int)std::floor(l / nb))) - 1;
          int J_MAX_ = (l - (nb * (int)std::floor(l / nb)) + 1) < nb
                           ? B * (l - (nb * (l / nb)) + 1)
                           : n - 1;
          int r_j = l < k ? l : l - k;
          int I_ = std::floor(r_j / nb) * B - 1 == -1 ? 0 : (r_j / nb) * B - 1;
          int I_MAX_ = ((r_j / nb) + 1) < nb ? ((r_j / nb) + 1) * B : n - 1;

          v.clear();
          v.resize((I_MAX_ - I_) * (J_MAX_ - J_));
          v = w.get_vector(I_, I_MAX_, J_, J_MAX_);
          MPI_Recv(&v[0], v.size(), MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);
          w.update_with_vector(I_, I_MAX_, J_, J_MAX_, v);
        }
        break;
      }

      iter++;
      // AQUI
      update_collumns(w, I, I_MAX, J, J_MAX, rank, n, nb);
      v.clear();
      v.resize((I_MAX - I) * (J_MAX - J));
      v = w.get_vector(I, I_MAX, J, J_MAX);
      MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD);
      u = w;
    }

  } else {
    for (;;) {
      if (rank < k) {
        a = red(I, I_MAX, J, J_MAX, rank, u, w);
        MPI_Recv(&b, 1, MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD, &status);
        a = b < a ? a : b;
        MPI_Send(&a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        v.clear();
        v.resize((I_MAX - I) * (J_MAX - J));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD,
                 &status);
        w.update_black_vector(I, I_MAX, J, J_MAX, v);
        //  MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (msg) {
          v.clear();
          v.resize((I_MAX - I) * (J_MAX - J));
          v = w.get_vector(I, I_MAX, J, J_MAX);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
          break;
        }
        update_collumns(w, I, I_MAX, J, J_MAX, rank, n, nb);
        v.clear();
        v.resize((I_MAX - I) * (J_MAX - J));
        v = w.get_vector(I, I_MAX, J, J_MAX);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD);

        u = w;
      } else {
        a = black(I, I_MAX, J, J_MAX, rank, u, w);

        MPI_Send(&a, 1, MPI_DOUBLE, rank - k, 0, MPI_COMM_WORLD);
        v = w.get_vector(I, I_MAX, J, J_MAX);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank - k, 0, MPI_COMM_WORLD);
        //  MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (msg) {
          break;
        }
        v.clear();
        v.resize((I_MAX - I) * (J_MAX - J));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank - k, 0, MPI_COMM_WORLD,
                 &status);
        // w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
        w.update_red_vector(I, I_MAX, J, J_MAX, v);
        w.update_sides(I, I_MAX, J, J_MAX, n, v);
        u = w;
      }
    }
  }
  return w;
}

void block_independent::poisson_gs(const int n) {
  int rank, total;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total);
  if (total != 2 * pow(std::floor(sqrt(total / 2)), 2)) {
    if (rank == 0) {
      std::cout << "Please run this program with " << std::floor(sqrt(total))
                << " MPI processes\n";
    }
    MPI_Finalize();
    exit(1);
  }

  Matrice<double> u(n, n);
  Matrice<double> w(n, n);
  Matrice<double> res(n, n);

  init_matrices(u, w, n);
  res = mpi_logic(u, w, res, n, rank, total);

  if (rank == 0) {
    std::cout << res << std::endl;
  }
  MPI_Finalize();
  // std::cout << iter << std::endl;
}
