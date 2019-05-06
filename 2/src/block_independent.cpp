#include "block_independent.h"

namespace {
double black(const int I, const int I_MAX, const int J, const int J_MAX,
             int rank, Matrice<double>& u, Matrice<double>& w) {
  int jstart;
  double max;
  double aux;

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
  double max;
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

Matrice<double>& mpi_logic(Matrice<double>& u, Matrice<double>& w,
                           Matrice<double>& res, int n, int rank, int NP) {
  double tol = pow((1.0 / n), 2);
  double diff = tol + 1;
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
        std::cout << "Recebi msg do rank " << l << "com o max " << b
                  << std::endl;
        a = a < b ? b : a;
      }

      diff = a;
      std::cout << "Max " << diff << std::endl;
      msg = diff <= tol;
      if (msg) {
        for (int l = 1; l < k; l++) {
          int J_ = B * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
                       ? 0
                       : B * (rank - (nb * (int)std::floor(rank / nb))) - 1;
          int J_MAX_ = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                           ? B * (rank - (nb * (rank / nb)) + 1)
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

      if (I != 0) {
        v = w.get_vector(I + 1, I + 2, J + 1, J_MAX - 1);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank - nb, 0, MPI_COMM_WORLD);
        v.clear();
        v.resize((J_MAX - J - 2));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank - nb, 0, MPI_COMM_WORLD,
                 &status);
        w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
      }
      if (I_MAX != n) {
        v = w.get_vector(I_MAX - 2, I_MAX - 1, J + 1, J_MAX - 1);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + nb, 0, MPI_COMM_WORLD);
        v.clear();
        v.resize((J_MAX - J - 2));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + nb, 0, MPI_COMM_WORLD,
                 &status);
        w.update_with_vector(I_MAX - 1, I_MAX, J + 1, J_MAX - 1, v);
      }
      if (J != 0) {
        v = w.get_vector(I + 1, I_MAX + 1, J + 1, J + 2);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        v.clear();
        v.resize((I_MAX - I - 2));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                 &status);
        w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
      }
      if (J_MAX != n) {
        v = w.get_vector(I + 1, I_MAX + 1, J_MAX - 2, J_MAX - 1);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        v.clear();
        v.resize((I_MAX - I - 2));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                 &status);
        w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
      }
      v.clear();
      v.resize((I_MAX - I) * (J_MAX - J));
      v = w.get_vector(I, I_MAX, J, J_MAX);
      MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD);
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
        std::cout << "Recebi do 0 "
                  << " " << msg << " " << rank << std::endl;
        if (msg) {
          v.clear();
          v.resize((I_MAX - I) * (J_MAX - J));
          v = w.get_vector(I, I_MAX, J, J_MAX);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
          break;
        }
        if (I != 0) {
          v = w.get_vector(I + 1, I + 2, J + 1, J_MAX - 1);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank - nb, 0, MPI_COMM_WORLD);
          v.clear();
          v.resize((J_MAX - J - 2));
          MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank - nb, 0, MPI_COMM_WORLD,
                   &status);
          w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
        }
        if (I_MAX != n) {
          v = w.get_vector(I_MAX - 2, I_MAX - 1, J + 1, J_MAX - 1);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + nb, 0, MPI_COMM_WORLD);
          v.clear();
          v.resize((J_MAX - J - 2));
          MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + nb, 0, MPI_COMM_WORLD,
                   &status);
          w.update_with_vector(I_MAX - 1, I_MAX, J + 1, J_MAX - 1, v);
        }
        if (J != 0) {
          v = w.get_vector(I + 1, I_MAX + 1, J + 1, J + 2);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
          v.clear();
          v.resize((I_MAX - I - 2));
          MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                   &status);
          w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
        }
        if (J_MAX != n) {
          v = w.get_vector(I + 1, I_MAX + 1, J_MAX - 2, J_MAX - 1);
          MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
          v.clear();
          v.resize((I_MAX - I - 2));
          MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                   &status);
          w.update_with_vector(I, I + 1, J + 1, J_MAX - 1, v);
        }
        v.clear();
        v.resize((I_MAX - I) * (J_MAX - J));
        v = w.get_vector(I, I_MAX, J, J_MAX);
        MPI_Send(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD);
      } else {
        a = black(I, I_MAX, J, J_MAX, rank, u, w);
        // mandar b para o processo 0
        //   send w to rank + k
        //   recieve w rank + k
        //   atualizar pretos
        //   se necessario iterar mais partilhar linhas
        // se nao mandar matrix para o 0
        //   exit
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
      }
    }
  }
  MPI_Finalize();
  return res;
}

void block_independent::poisson_gs(int argc, char* argv[], const int n) {
  int rank, total;
  std::cout << "in mpi_logic" << std::endl;
  MPI_Init(&argc, &argv);
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
  mpi_logic(u, w, res, n, rank, total);

  std::cout << res << std::endl;
  // std::cout << iter << std::endl;
}
