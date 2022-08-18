/*

Pade approximation

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _PADE_H
#define _PADE_H

#include <complex>

//
// Analytic continuation by the Pade approximation
//
// u({z_i}) : the original discrete function.
// a : coefficient set of the approximated function 'u_{pade}(z)'
// N : length of arrays 'z', 'u' and 'a'
//

// initialize 'a'
void pade_init(std::complex<double> *z, std::complex<double> *u, std::complex<double> *a, int N);

// calculate approximated value 'u_{pade}(w)' at 'w'
std::complex<double> pade_calc(std::complex<double> *z, std::complex<double> *a, std::complex<double> w, int N);

#endif // _PADE_H
