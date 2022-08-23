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

#include "mt.h"
#include "fft.h"
#include "gtau.h"


#ifndef _IMAG
#define _IMAG
const std::complex<double> IMAG(0, 1.0);
#endif // _IMAG


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



#endif // _CT_QMC_SHARE_H
