#ifndef _VECTOR_TYPE_H
#define _VECTOR_TYPE_H

#include <vector>
#include <complex>

using vec_d = std::vector<double>;
using vec_c = std::vector<std::complex<double> >;
using vec_vec_d = std::vector<std::vector<double> >;
using vec_vec_c = std::vector<std::vector<std::complex<double> > >;
using vec_vec_vec_d = std::vector<std::vector<std::vector<double> > >;
using vec_vec_vec_c = std::vector<std::vector<std::vector<std::complex<double> > > >;
using vec_vec_vec_vec_d = std::vector<std::vector<std::vector<std::vector<double> > > >;
using vec_vec_vec_vec_c = std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > >;


template<class T>
void zeros(std::vector<T> &x)
{
    std::fill(x.begin(), x.end(), 0);
}

template<class T>
void zeros(std::vector<std::vector<T> > &x)
{
    for (auto& x1 : x) {
        std::fill(x1.begin(), x1.end(), 0);
    }
}

template<class T>
void zeros(std::vector<std::vector<std::vector<T> > > &x)
{
    for (auto& x1 : x) {
        for (auto& x2 : x1) {
            std::fill(x2.begin(), x2.end(), 0);
        }
    }
}

template<class T>
void zeros(std::vector<std::vector<std::vector<std::vector<T> > > > &x)
{
    for (auto& x1 : x) {
        for (auto& x2 : x1) {
            for (auto& x3 : x2) {
                std::fill(x3.begin(), x3.end(), 0);
            }
        }
    }
}

template<class T>
void resize(std::vector<std::vector<T> > &x, std::size_t n1, std::size_t n2)
{
    x.resize(n1);
    for (auto& x1 : x) {
        x1.resize(n2);
    }
    // for (size_t i=0; i<n1; i++) {
    //     x[i].resize(n2);
    // }
}

template<class T>
void resize(std::vector<std::vector<std::vector<T> > > &x, std::size_t n1, std::size_t n2, std::size_t n3)
{
    x.resize(n1);
    for (auto& x1 : x) {
        x1.resize(n2);
        for (auto& x2 : x1) {
            x2.resize(n3);
        }
    }
    // for (size_t i=0; i<n1; i++) {
    //     x[i].resize(n2);
    //     for (size_t j=0; j<n2; j++) {
    //         x[i][j].resize(n3);
    //     }
    // }
}

template<class T>
void resize(std::vector<std::vector<std::vector<std::vector<T> > > > &x, std::size_t n1, std::size_t n2, std::size_t n3, std::size_t n4)
{
    x.resize(n1);
    for (auto& x1 : x) {
        x1.resize(n2);
        for (auto& x2 : x1) {
            x2.resize(n3);
            for (auto& x3 : x2) {
                x3.resize(n4);
            }
        }
    }
}


template<class T>
bool check_size(const std::vector<T> &x, std::size_t n1)
{
    return x.size() == n1;
}

template<class T>
bool check_size(const std::vector<std::vector<T> > &x, std::size_t n1, std::size_t n2)
{
    bool r = (x.size() == n1);
    for (auto& x1 : x) {
        r &= (x1.size() == n2);
    }
    return r;
}

template<class T>
bool check_size(const std::vector<std::vector<std::vector<T> > > &x, std::size_t n1, std::size_t n2, std::size_t n3)
{
    bool r = (x.size() == n1);
    for (auto& x1 : x) {
        r &= (x1.size() == n2);
        for (auto& x2 : x1) {
            r &= (x2.size() == n3);
        }
    }
    return r;
}

template<class T>
bool check_size(const std::vector<std::vector<std::vector<std::vector<T> > > > &x, std::size_t n1, std::size_t n2, std::size_t n3, std::size_t n4)
{
    bool r = (x.size() == n1);
    for (auto& x1 : x) {
        r &= (x1.size() == n2);
        for (auto& x2 : x1) {
            r &= (x2.size() == n3);
            for (auto& x3 : x2) {
                r &= (x3.size() == n4);
            }
        }
    }
    return r;
}

#endif // _VECTOR_TYPE_H
