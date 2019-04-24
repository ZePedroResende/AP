#ifndef MATRICE_H
#define MATRICE_H
#include <memory>

template <class T>
class Matrice {
    int m_width;
    int m_height;
    std::unique_ptr<T[]> array;

   public:
    Matrice(int width, int height)
        : m_width(width), m_height(height), array(std::make_unique<T[]>(width * height)) {}
    T& operator()(int x, int y) { return array[x + m_width * y]; };
};
#endif  // MATRICE_H
