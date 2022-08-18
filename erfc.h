#ifndef _ERFC_H
#define _ERFC_H

#include <complex>

//
// return W(z)
// W(z) = -i sqrt(pi) w(z)
// w(z) = exp(-z^2) erfc(-iz)
//
std::complex<double> erfc_w(std::complex<double> z);

#endif // _ERFC_H
