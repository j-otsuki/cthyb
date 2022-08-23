/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _GTAU_H
#define _GTAU_H

// #include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <complex>

struct green_func_0{
	int N;  // N_TAU
	gsl_interp_accel *acc_p;
	gsl_interp_accel *acc_m;
	gsl_spline *spline_p;  // G0(tau)  [0:beta]
	gsl_spline *spline_m;  // G0(-tau)  [0:beta]
};

struct boson_func_0{
	int N;  // N_TAU
	gsl_interp_accel *acc_p;
	gsl_interp_accel *acc_m;
	gsl_spline *spline_p;
	gsl_spline *spline_m;
};


// interporation algorithm
#define INTERP gsl_interp_linear
#define INTERP_CHR "gsl_interp_linear"

// #define INTERP gsl_interp_cspline
// #define INTERP_CHR "gsl_itnerp_cspline"

// #define INTERP gsl_interp_akima
// #define INTERP_CHR "gsl_interp_akima"



// void G0_alloc(struct green_func_0 &G0);
void G0_alloc(struct green_func_0 &G0, int N);
void G0_free(struct green_func_0 &G0);

// void G0_init_integ(struct green_func_0 &G0, double beta, double D);
// void G0_init_integ(struct green_func_0 &G0, double beta, double D, double V_sqr);

void G0_init(struct green_func_0 &G0, double *G_tau, double beta);

// a = iw * G(iw), w->inf
void G0_init_fft(struct green_func_0 &G0, std::complex<double> *G_omega, double beta, double a);
// a=1
void G0_init_fft(struct green_func_0 &G0, std::complex<double> *G_omega, double beta);

// return G0(tau2-tau1)
// if tau2-tau1==0, return G0(+0)
double G0_calc_interp(struct green_func_0 &G0, double tau2, double tau1);
// if tau2-tau1==0, return G0(-0)
double G0_calc_interp_0m(struct green_func_0 &G0, double tau2, double tau1);

double G0_calc_integ(double tau, double beta, double D);



void D0_alloc(struct boson_func_0 &D0, int N);

void D0_free(struct boson_func_0 &D0);

// a = iw * D(iw), w->inf
void D0_init_fft(struct boson_func_0 &D0, std::complex<double> *D_omega, double beta, double a);
// a = 0
void D0_init_fft(struct boson_func_0 &D0, std::complex<double> *D_omega, double beta);

// return D0(tau2-tau1)
// if tau2-tau1==0, return D0(+0)
double D0_calc_interp(struct boson_func_0 &D0, double tau2, double tau1);
// if tau2-tau1==0, return D0(-0)
double D0_calc_interp_0m(struct boson_func_0 &D0, double tau2, double tau1);

void K0_init_sum(struct boson_func_0 &K0, std::complex<double> *D_omega, double beta);
void K1_init_sum(struct boson_func_0 &K1, std::complex<double> *D_omega, double beta);


#endif // _GTAU_H
