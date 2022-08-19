/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _CT_QMC_SHARE_H
#define _CT_QMC_SHARE_H

#include <stdio.h>
#include <math.h>
#include <complex>
#include <ctime>
#include <string.h>
#include <vector>
// #include <gsl/gsl_integration.h>
// #include <gsl/gsl_spline.h>
// #include <gsl/gsl_sf_erf.h>
// #include "green_func_0.h"
#include "mt.h"
#include "fft.h"
// #include "erfc.h"
// #include "pade.h"
#include "gtau.h"

// using namespace std;

// #define N_K 1024  // maximum number of k
// #define N_TAU 1024  // discretization of beta
#define DATA_DIR "output/"

// #define DOS 0
// 0: constant
// 1: semicircle (Bethe lattice)
// 2: Gaussian (hyper-cubic lattice)
// 3: d=2 tight-binding band (N_L=256, 512)
// 4: d=3 tight-binding band (N_L=32, 64)
// 5: d=4 tight-binding band (N_L=16, 32)
// 6: d=5 tight-binding band (N_L=8, 16)
// -1: single bath

// #define N_L 32  // # of k-points between k=0 and PI

const double ACCURACY_FILLING = 1e-8;
const int FILLING_ITER_MAX = 100;

// #define EQUAL_TIME_G 0
// 0 for antiferro-exchange
// 1 for ferro-exchange interaction (only for N_F=1 and 2)
// 1 for U-expansion



#ifndef _IMAG
#define _IMAG
const std::complex<double> IMAG(0, 1.0);
#endif // _IMAG

// #define MAX(a,b)  ((a) > (b) ? (a) : (b))
// #define MIN(a,b)  ((a) > (b) ? (b) : (a))


#define _LAPACK
#ifdef _LAPACK
extern "C"{
	void dgetrf_(long*, long*, double*, long*, long*, long*);
	void dgetri_(long*, double*, long*, long*, double*, long*, long*);
	void zgetrf_(long*, long*, std::complex<double>*, long*, long*, long*);
	void zgetri_(long*, std::complex<double>*, long*, long*, std::complex<double>*, long*, long*);
}
#endif // _LAPACK



// // int [0:n)
// inline int rand_int(int n);
// // [0:beta)
// inline double rand_tau_beta(double beta);
// // (0:l_max)
// inline double rand_tau_l(double l_max);
// // int [0:n)

// int [0:n)
inline int rand_int(int n)
{
	return( genrand_int32() % n );
}

// [0:beta)
inline double rand_tau_beta(double beta)
{
// 	return( genrand_real2()*beta );
	return( genrand_res53()*beta );
}

// (0:l_max)
inline double rand_tau_l(double l_max)
{
	return( genrand_real3()*l_max );
}

void random_permutation(int *array, int n);


double sgn(double);
double pow_m1(int);


// 1 / (e^x + 1) = ( 1-tanh(x/2) ) / 2
double distrib_fermi(double x);

// e^y / (e^x + 1)
double distrib_fermi_boltz(double x, double y);


int tau_order(double *, int, double);
int tau_order(std::vector<double>&, double);
int tau_position(double *, int, double);


//
// 0: do not update,  1: update
//
// reject if prob<0
int metropolis(double prob, unsigned long &n_accept, unsigned long &n_reject);
// input fabs(prob)
int metropolis_abs(double prob, unsigned long &n_accept, unsigned long &n_reject);


// return: number of data stored in params[]
int read_data(FILE *fp, double *params);

void print_time(clock_t time1, clock_t time2);
void print_time(clock_t time1, clock_t time2, char *filename);
void print_time_mpi(double time_trans);
void print_time_mpi(clock_t time1, clock_t time2, double time_trans, char *filename);
void sprint_time(char *str, clock_t time1);
void sprint_time_mpi(char *str, clock_t time1, double time_trans);
void sprint_time_mpi(char *str, clock_t time1, clock_t time2, double time_trans);


void mesh_init(int i_mesh, double *y, double x_min, double x_max, int N);
void mesh_init(int i_mesh, double *y, double x_min, double x_max, int N, int my_rank);
void mesh_linear(double *y, double x_min, double x_max, int N);
void mesh_log_linear(double *y, double x_min, double x_max, int N);

//
// non-linear mesh
// (n_tau1+1) points in the range [0:beta/2]
// shortest interval: beta/2/n_tau2
void tau_mesh_nonlinear_boson(double *tau, int n_tau1, int n_tau2, double beta);
//
// (n_tau1+1) points in the range [0:beta]
// shortest interval: beta/n_tau2
void tau_mesh_nonlinear_fermion(double *tau, int n_tau1, int n_tau2, double beta);


// struct cond_op{
// 	int k;  // number of segments
// 	double tau1[N_K];  // for f-annihilation (c-creation) operator
// 	double tau2[N_K];  // for f-creation (c-annihilation) operator
// 	double mat_M[N_K][N_K];  // [tau1][tau2]
// 	int flag;  // 0: tau1[i] < tau2[i],  1: tau1[i] > tau2[i]
// };


void phase_init(struct phase &PHASE);
std::complex<double> phase_interp(struct phase &PHASE, double x);
void phase_free(struct phase &PHASE);


//
// for test
//
void print_mat_M(struct cond_op &F);

void test_calc_mat_Delta(struct cond_op &F, struct green_func_0 &G0);
void test_calc_mat_M(struct cond_op &F, struct green_func_0 &G0);
void test_product_Delta_M(struct cond_op &F, struct green_func_0 &G0);


#endif // _CT_QMC_SHARE_H
