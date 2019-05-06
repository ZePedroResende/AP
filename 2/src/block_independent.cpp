#include "block_independent.h"

namespace {
double black(const int max_i, const int max_j, int rank, const int k,
             Matrice<double>& u, Matrice<double>& w) {
  std::cout << rank << " entrei no black" << std::endl;
  int jstart;
  double max;
  double aux;
  rank -= k;
  int nb = max_i / k;
  int b = 2 * nb;
  int I = b * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
              ? 0
              : b * (rank - (nb * (int)std::floor(rank / nb))) - 1;
  int I_MAX = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                  ? b * (rank - (nb * (rank / nb)) + 1)
                  : max_i - 1;
  int J = std::floor(rank / nb) * b - 1 == -1 ? 0 : (rank / nb) * b - 1;
  int J_MAX = ((rank + 1 / nb)) < nb ? ((rank + 1 / nb)) * b : max_j - 1;
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
  std::cout << rank + k << " sai do black" << std::endl;
  return max;
}

double red(const int max_i, const int max_j, const int rank, const int k,
           Matrice<double>& u, Matrice<double>& w) {
  std::cout << rank << " entrei no red" << std::endl;
  int jstart;
  double max;
  double aux;
  int nb = max_i / k;
  int b = 2 * nb;
  int I = b * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
              ? 0
              : b * (rank - (nb * (int)std::floor(rank / nb))) - 1;
  int I_MAX = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                  ? b * (rank - (nb * (rank / nb)) + 1)
                  : max_i - 1;
  int J = std::floor(rank / nb) * b - 1 == -1 ? 0 : (rank / nb) * b - 1;
  int J_MAX = ((rank / nb) + 1) < nb ? ((rank / nb) + 1) * b : max_j - 1;
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
  std::cout << rank << " sai do red" << std::endl;
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
  double a;
  int iter = 0;
  // mpi
  int P = NP / 2;
  int X = std::floor(sqrt(P));
  int k = pow(X, 2);
  int msg;
  MPI_Status status;

  std::vector<double> v;
  int nb = std::floor(sqrt(k));
  int b = std::floor(n / sqrt(k));
  int I = b * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
              ? 0
              : b * (rank - (nb * (int)std::floor(rank / nb))) - 1;
  int I_MAX = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                  ? b * (rank - (nb * (rank / nb)) + 1)
                  : n - 1;
  int rank_j = rank < k ? rank : rank - k;
  int J = std::floor(rank_j / nb) * b - 1 == -1 ? 0 : (rank_j / nb) * b - 1;
  int J_MAX = ((rank_j / nb) + 1) < nb ? ((rank_j / nb) + 1) * b : n - 1;
  std::cout << "rank=" << rank << " nb=" << nb << " b=" << b << " I = " << I
            << " I_MAX=" << I_MAX << " J = " << J << " J_MAX=" << J_MAX
            << std::endl;

  if (rank == 0) {
    for (;;) {
      a = red(n - 1, n - 1, rank, k, u, w);
      MPI_Recv(&msg, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
      std::cout << "Recebi msg do meu black " << rank << " !" << std::endl;
      a = msg < a ? a : msg;

      v.clear();
      v.resize((I_MAX - I) * (J_MAX - J));
      MPI_Recv(&v[0], v.size(), MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
      w.update_black_vector(I, I_MAX, J, J_MAX, v);

      for (int l = 1; l < k; l++) {
        b = MPI_Recv(&msg, 1, MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);
        a = a < b ? b : a;
      }

      diff = a;
      msg = diff <= tol;
      MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (msg) {
        for (int l = 1; l < k; l++) {
          int I_ = b * (rank - (nb * (int)std::floor(rank / nb))) - 1 == -1
                       ? 0
                       : b * (rank - (nb * (int)std::floor(rank / nb))) - 1;
          int I_MAX_ = (rank - (nb * (int)std::floor(rank / nb)) + 1) < nb
                           ? b * (rank - (nb * (rank / nb)) + 1)
                           : n - 1;
          int r_j = l < k ? l : l - k;
          int J_ = std::floor(r_j / nb) * b - 1 == -1 ? 0 : (r_j / nb) * b - 1;
          int J_MAX_ = ((r_j / nb) + 1) < nb ? ((r_j / nb) + 1) * b : n - 1;

          v.clear();
          v.resize((I_MAX_ - I_) * (J_MAX_ - J_));
          v = w.get_vector(I_, I_MAX_, J_, J_MAX_);
          b = MPI_Recv(&v[0], v.size(), MPI_DOUBLE, l, 0, MPI_COMM_WORLD,
                       &status);
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
        a = red(n - 1, n - 1, rank, k, u, w);
        MPI_Recv(&msg, 1, MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD, &status);
        std::cout << "Recebi msg do meu black " << rank << " !" << std::endl;
        a = msg < a ? a : msg;
        MPI_Send(&a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        v.clear();
        v.resize((I_MAX - I) * (J_MAX - J));
        MPI_Recv(&v[0], v.size(), MPI_DOUBLE, rank + k, 0, MPI_COMM_WORLD,
                 &status);
        w.update_black_vector(I, I_MAX, J, J_MAX, v);
        MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        std::cout << msg << rank << std::endl;
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
        a = black(n - 1, n - 1, rank, k, u, w);
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
        MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
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
