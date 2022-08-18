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


template<class T>
void zeros(std::vector<T> &x)
{
    for (const auto& ele : x) {
        ele = 0;
    }
}

template<class T>
void zeros(std::vector<std::vector<T> > &x)
{
    for (const auto& row : x) {
        for (const auto& ele : row) {
            ele = 0;
        }
    }
}

template<class T>
void resize(std::vector<std::vector<T> > &x, int n1, int n2)
{
    x.resize(n1);
    // for (const auto& x1 : x) {
    for (auto& x1 : x) {
        x1.resize(n2);
    }
}

template<class T>
void resize(std::vector<std::vector<T> > &x, int n1, int n2, int n3)
{
    x.resize(n1);
    for (auto& x1 : x) {
        x1.resize(n2);
        for (auto& x2 : x1) {
            x2.resize(n3);
        }
    }
}

#endif // _VECTOR_TYPE_H
