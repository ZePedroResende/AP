#ifndef MATRICE_H
#define MATRICE_H
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

template <class T>
class Matrice {
 public:
  int m_height;
  int m_width;
  std::vector<T> array;

  constexpr Matrice(const int height, const int width)
      : m_height(height),
        m_width(width),
        array(std::vector<T>(width * height)) {}

  constexpr Matrice(const Matrice<T>& m)
      : m_height(m.m_height),
        m_width(m.m_width),
        array(std::vector<T>(m.array)) {}

  T max() { return *std::max_element(array.begin(), array.end()); }

  constexpr T& operator()(const int i, const int j) {
    return array[j + (m_width * i)];
  };

  Matrice<T> operator-(Matrice<T> b) {
    Matrice<T> a = *this;
    Matrice<T> res(a.m_height, a.m_width);

    for (int i = 0; i < res.m_height; i++) {
      for (int j = 0; j < res.m_width; j++) {
        res(i, j) = a(i, j) - b(i, j);
      }
    }

    return res;
  }

  friend std::ostream& operator<<(std::ostream& os, Matrice<T>& m) {
    for (int i = 0; i < m.m_height; i++) {
      for (int j = 0; j < m.m_width; j++) {
        os << m(i, j) << " ";
      }
      os << std::endl;
    }
    return os;
  }

  std::vector<T> get_vector(int I, int I_MAX, int J, int J_MAX) {
    std::vector<T> res;

    Matrice<T>& a = *this;
    for (int i = I; i < I_MAX; i++) {
      for (int j = J; j < J_MAX; j++) {
        res.push_back(a(i, j));
      }
    }

    return res;
  }

  bool update_with_vector(int I, int I_MAX, int J, int J_MAX,
                          std::vector<T>& in) {
    Matrice<T>& a = *this;
    int it = 0;
    for (int i = I; i < I_MAX; i++) {
      for (int j = J; j < J_MAX; j++) {
        if (I == 2 && j == 3) {
          /*
          std::cout << i << " " << j << " " << it << " " << a(i, j) << " "
                    << in[it] << "\n\n";
                    */
        }
        // std::cout << I << I_MAX << J << J_MAX << it << std::endl;
        a(i, j) = in[it];
        it++;
      }
    }

    return in.size() == it;
  }

  void update_black_vector(int I, int I_MAX, int J, int J_MAX,
                           std::vector<double>& in) {
    int jstart;
    Matrice<T>& w = *this;
    for (int i = (I + 1), i_ = 0; i < I_MAX; i++, i_++) {
      if (i % 2 == 1)
        jstart = 2;  // odd row
      else
        jstart = 1;  // even row
      for (int j = jstart + J, j_ = jstart; j < J_MAX; j += 2, j_ += 2) {
        w(i, j) = in[j_ + ((J_MAX - J + 1) * i_)];
      }
    }
  }

  void update_red_vector(int I, int I_MAX, int J, int J_MAX,
                         std::vector<double>& in) {
    int jstart;
    Matrice<T>& w = *this;
    for (int i = (I + 1), i_ = 0; i < I_MAX; i++, i_++) {
      if (i % 2 == 1)
        jstart = 1;  // odd row
      else
        jstart = 2;  // even row
      for (int j = jstart + J, j_ = jstart; j < J_MAX; j += 2, j_ += 2) {
        w(i, j) = in[j_ + ((J_MAX - J + 1) * i_)];
      }
    }
  }

  void update_sides(int I, int I_MAX, int J, int J_MAX, int n,
                    std::vector<double>& v) {
    Matrice<T>& w = *this;

    if (I != 0) {
      w.update_with_vector(I, I + 1, J + 1, J_MAX, v);
    }
    if (I_MAX != n) {
      w.update_with_vector(I_MAX - 1, I_MAX, J + 1, J_MAX, v);
    }
    if (J != 0) {
      w.update_with_vector(I + 1, I_MAX, J, J + 1, v);
    }
    if (J_MAX != n) {
      w.update_with_vector(I + 1, I_MAX, J_MAX - 1, J_MAX, v);
    }
  }
};
#endif  // MATRICE_H
