/*
Contiguous multidimensional array

Example:
    Array2D<double> a(n1, n2);
    double* b = a.data();

    Then access by
        a(i, j)
    and
        b[i*n2+j]
    are equivalent.
*/

#ifndef _ARRAY_HPP
#define _ARRAY_HPP

#include <vector>


template <typename T>
class Array2D{
public:
    Array2D() {};
    Array2D(std::size_t n1, std::size_t n2) : n1(n1), n2(n2), array(n1*n2) {};

    // a(i, j)
    T& operator ()(int i, int j){
        return array[i * n2 + j];
    }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return array.size(); }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array2D& array) { array.zeros(); }

private:
    std::size_t n1, n2;
    std::vector<T> array;
};


template <typename T>
class Array3D{
public:
    Array3D() {};
    Array3D(std::size_t n1, std::size_t n2, std::size_t n3) : n1(n1), n2(n2), n3(n3), array(n1*n2*n3) {};

    // a(i, j, k)
    T& operator ()(int i, int j, int k){
        return array[(i * n2 + j) * n3 + k];
    }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return array.size(); }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array3D& array) { array.zeros(); }

private:
    std::size_t n1, n2, n3;
    std::vector<T> array;
};


template <typename T>
class Array4D{
public:
    Array4D() {};
    Array4D(std::size_t n1, std::size_t n2, std::size_t n3, std::size_t n4) : n1(n1), n2(n2), n3(n3), n4(n4), array(n1*n2*n3*n4) {};

    // a(i, j, k, l)
    T& operator ()(int i, int j, int k, int l){
        return array[((i * n2 + j) * n3 + k) * n4 + l];
    }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return array.size(); }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array4D& array) { array.zeros(); }

private:
    std::size_t n1, n2, n3, n4;
    std::vector<T> array;
};


template <typename T>
class Array5D{
public:
    Array5D() {};
    Array5D(std::size_t n1, std::size_t n2, std::size_t n3, std::size_t n4, std::size_t n5) : n1(n1), n2(n2), n3(n3), n4(n4), n5(n5), array(n1*n2*n3*n4*n5) {};

    // a(i, j, k, l, m)
    T& operator ()(int i, int j, int k, int l, int m){
        return array[(((i * n2 + j) * n3 + k) * n4 + l) * n5 + m];
    }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return array.size(); }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array5D& array) { array.zeros(); }

private:
    std::size_t n1, n2, n3, n4, n5;
    std::vector<T> array;
};


#endif // _ARRAY_HPP
