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

    bool update_with_vector_agregate(int I, int I_MAX, int J, int J_MAX,
        std::vector<T>& in) {
      Matrice<T>& a = *this;
      int it = 0;
      for (int i = I; i < I_MAX; i++) {
        for (int j = J; j < J_MAX; j++) {
          if(in[it])
            a(i, j) = in[it];
            it++;
        }
      }

      return in.size() == it;
    }
    bool update_with_vector(int I, int I_MAX, int J, int J_MAX,
        std::vector<T>& in) {
      Matrice<T>& a = *this;
      int it = 0;
      for (int i = I; i < I_MAX; i++) {
        for (int j = J; j < J_MAX; j++) {
          a(i, j) = in[it];
          it++;
        }
      }

      return in.size() == it;
    }
    void update_with_v(int I, int I_MAX, int J, int J_MAX,
        std::vector<T>& in, int I_, int I_MAX_, int J_, int J_MAX_) {
      Matrice<T>& a = *this;

      /*
         std::cout <<"\n" ;
         for (int i = I, i_ = I_; i < I_MAX; i++, i_++) {
         for (int j = J, j_ = J_; j < J_MAX; j++, j_++) {
         std::cout <<" " << in[j_ + ((J_MAX_ - J_ + 1) * i_)];
         }
         std::cout <<"\n" ;
         }
         std::cout <<"\n" ;
         */

      for (int i = I, i_ = I_; i < I_MAX; i++, i_++) {
        for (int j = J, j_ = J_; j < J_MAX; j++, j_++) {
          a(i, j) = in[j_ + ((J_MAX_ - J_ + 1) * i_)];
        }
      }

    }

    void update_black_vector(int I, int I_MAX, int J, int J_MAX,
        std::vector<double>& in) {
      int jstart;
      Matrice<T>& w = *this;
      for (int i = (I + 1), i_ = 1; i < I_MAX; i++, i_++) {
        if (i % 2 == 1)
          jstart = 2;  // odd row
        else
          jstart = 1;  // even row
        for (int j = jstart + J, j_ = jstart; j < J_MAX; j += 2, j_ += 2) {
          w(i, j) = in[j_ + ((J_MAX - J ) * i_)];
        }
      }
    }

    void update_red_vector(int I, int I_MAX, int J, int J_MAX,
        std::vector<double>& in) {
      int jstart;
      Matrice<T>& w = *this;
      for (int i = (I + 1), i_ = 1; i < I_MAX; i++, i_++) {
        if (i % 2 == 1)
          jstart = 1;  // odd row
        else
          jstart = 2;  // even row
        for (int j = jstart + J, j_ = jstart; j < J_MAX; j += 2, j_ += 2) {
          w(i, j) = in[j_ + ((J_MAX - J ) * i_)];
        }
      }
    }

    void update_sides(int I, int I_MAX, int J, int J_MAX, int n,
        std::vector<double>& v) {
      Matrice<T>& w = *this;

      int i, i_max, j, j_max;
      int i_, i_max_, j_, j_max_;


      i_ = I == 0 ? I : I + 1;
      i_max_ = I_MAX == n - 1 ? I_MAX : I_MAX - 1;
      j_ = J == 0 ? J : J + 1;
      j_max_ = J_MAX == n - 1 ? J_MAX : J_MAX - 1;


      i = I == 0 ? 0 :  1;
      i_max = I_MAX == n - 1 ? i_max_ - i_ + 1 : i_max_ - i_ ;
      j = J == 0 ? 0 :  1;
      j_max = J_MAX == n - 1 ? j_max_ - j_ + 1 : j_max_ - j_;


      if (I != 0) {
        w.update_with_v(I, I + 1, j_, j_max_ + 1, v, 0,1,j,j_max_ - j + 1);
      }
      if (I_MAX != n -1) {
        w.update_with_v(I_MAX , I_MAX + 1, j_, j_max_ + 1, v, I_MAX - I, I_MAX - I + 1,j, j_max + 1);
      }
      if (J != 0) {
        w.update_with_v(i_, i_max_ + 1, J, J + 1, v, i,i_max +1 , 0,1 );
      }
      if (J_MAX != n - 1) {
        w.update_with_v(i_, i_max_ + 1, J_MAX , J_MAX + 1 , v,i,i_max + 1, J_MAX - J, J_MAX -J +1 );
      }

    }

    bool replace(int I, int I_MAX, int J, int J_MAX,
        std::vector<T>& in) {
      Matrice<T>& a = *this;
      int it = 0;
      for (int i = I; i < I_MAX; i++) {
        for (int j = J; j < J_MAX; j++) {
          a(i, j) = in[it];
          it++;
        }
      }

      return in.size() == it;
    }

};
#endif  // MATRICE_H
