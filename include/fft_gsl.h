/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _FFT_GSL_H
#define _FFT_GSL_H

#include <complex>


// FFT for fermionic G

// double G_tau[N]
// complex<double> G_omega[N/2]

//
// a = G(-0) - G(+0)
//

// a=1
void fft_fermion_radix2_tau2omega(double *G_tau, std::complex<double> *G_omega, double beta, int N);
// a: arbitrary
void fft_fermion_radix2_tau2omega(double *G_tau, std::complex<double> *G_omega, double beta, int N, double a);

// a=1
void fft_fermion_radix2_omega2tau(double *G_tau, std::complex<double> *G_omega, double beta, int N);
// a: arbitrary
void fft_fermion_radix2_omega2tau(double *G_tau, std::complex<double> *G_omega, double beta, int N, double a);


// FFT for bosonic G

// double G_tau[N]
// complex<double> G_omega[N/2+1]

//
// a = G(-0) - G(+0)
//

// a=0
void fft_boson_radix2_tau2omega(double *G_tau, std::complex<double> *G_omega, double beta, int N);
// a: arbitrary
void fft_boson_radix2_tau2omega(double *G_tau, std::complex<double> *G_omega, double beta, int N, double a);

// a=0
void fft_boson_radix2_omega2tau(double *G_tau, std::complex<double> *G_omega, double beta, int N);
// a: arbitrary
void fft_boson_radix2_omega2tau(double *G_tau, std::complex<double> *G_omega, double beta, int N, double a);

#endif // _FFT_GSL_H
