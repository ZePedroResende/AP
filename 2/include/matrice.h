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

    Matrice(int height, int width)
        : m_height(height), m_width(width), array(std::vector<T>(width * height)) {}

    Matrice(const Matrice<T>& m)
        : m_height(m.m_height), m_width(m.m_width), array(std::vector<T>(m.array)) {}

    T max() { return *std::max_element(array.begin(), array.end()); }

    T& operator()(int x, int y) { return array[x + m_width * y]; };

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
};
#endif  // MATRICE_H
