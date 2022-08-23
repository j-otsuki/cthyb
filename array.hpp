/*
Contiguous multidimensional array

Example:
    Array2D<double> a(n1, n2);
    double* b = a.data();

    Then access by
        a[i][j]
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
    Array2D(std::size_t n1, std::size_t n2) : n1(n1), n2(n2), array(n1*n2), p_array(n1) {
        for(int i=0; i<n1; i++){
            p_array[i] = &array[i*n2];
        }
    }

    // a[i][j]
    T* operator [](int n) { return p_array[n]; }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return n1*n2; }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array2D& array) { array.zeros(); }

private:
    std::size_t n1, n2;
    std::vector<T> array;
    std::vector<T*> p_array;
};


template <typename T>
class Array3D{
public:
    Array3D() {};
    Array3D(std::size_t n1, std::size_t n2, std::size_t n3) : n1(n1), n2(n2), n3(n3), array(n1*n2*n3), pp_array(n1) {
        for(int i=0; i<n1; i++){
            std::vector<T*> p_array(n2);
            for(int j=0; j<n2; j++){
                p_array[j] = &array[(i*n2+j)*n3];
            }
            pp_array[i] = p_array;
        }
    }

    // a[i][j][k]
    std::vector<T*> operator [](int n) { return pp_array[n]; }

    // return pointer to the head address
    T* data() { return array.data(); }

    std::size_t size() { return n1*n2*n3; }

    void zeros() {
        for (auto &a : array)  a = 0;
    }
    friend void zeros(Array3D& array) { array.zeros(); }

private:
    std::size_t n1, n2, n3;
    std::vector<T> array;
    std::vector<std::vector<T*> > pp_array;
};


#endif // _ARRAY_HPP
