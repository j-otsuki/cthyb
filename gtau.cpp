/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "gtau.h"
// #include "green_func_0.h"
// #include <gsl/gsl_roots.h>
// #include <gsl/gsl_errno.h>
#include <complex>
#include "fft.h"

// #include "stdio.h"
// #define DATA_DIR "data/"

// #define EQUAL_TIME_G 0
// 0 for antiferro-exchange
// 1 for ferro-exchange interaction (only for N_F=1 and 2)
// 1 for U-expansion


using namespace std;

void G0_alloc(struct green_func_0 &G0, int N)
{
	G0.N = N;
	G0.acc_p= gsl_interp_accel_alloc ();
	G0.acc_m= gsl_interp_accel_alloc ();
	G0.spline_p = gsl_spline_alloc (INTERP, G0.N+1);  // Delta(tau)  [0:beta]
	G0.spline_m = gsl_spline_alloc (INTERP, G0.N+1);  // Delta(-tau)  [0:beta]
}

void G0_free(struct green_func_0 &G0)
{
	gsl_spline_free (G0.spline_p);
	gsl_spline_free (G0.spline_m);
	gsl_interp_accel_free (G0.acc_p);
	gsl_interp_accel_free (G0.acc_m);
}

// void G0_init_integ(struct green_func_0 &G0, double beta, double D)
// {
// 	double G0_calc_integ(double, double, double);
// 	double tau[G0.N+1], G0_p[G0.N+1], G0_m[G0.N+1];

// 	for(int i=0; i<=G0.N; i++){
// 		tau[i] = (double)i * beta / (double)G0.N;
// 		G0_p[i] = G0_calc_integ(tau[i], beta, D);
// 		G0_m[G0.N-i] = -G0_p[i];
// 	}

// 	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
// 	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);
// }
// void G0_init_integ(struct green_func_0 &G0, double beta, double D, double V_sqr)
// {
// 	double G0_calc_integ(double, double, double);
// 	double tau[G0.N+1], G0_p[G0.N+1], G0_m[G0.N+1];

// 	for(int i=0; i<=G0.N; i++){
// 		tau[i] = (double)i * beta / (double)G0.N;
// 		G0_p[i] = V_sqr * G0_calc_integ(tau[i], beta, D);
// 		G0_m[G0.N-i] = -G0_p[i];
// 	}

// 	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
// 	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);
// }

void G0_init(struct green_func_0 &G0, double *G_tau, double beta)
{
	double *tau = new double[G0.N+1];
	double *G0_p = new double[G0.N+1];
	double *G0_m = new double[G0.N+1];

	for(int i=0; i<=G0.N; i++){
		tau[i] = (double)i * beta / (double)G0.N;
		G0_p[i] = G_tau[i];
		G0_m[G0.N-i] = -G0_p[i];
	}

	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);

	delete [] tau;
	delete [] G0_p;
	delete [] G0_m;
}

// a = iw * G(iw), w->inf
void G0_init_fft(struct green_func_0 &G0, complex<double> *G_omega, double beta, double a)
{
// 	double tau[G0.N+1], G0_p[G0.N+1], G0_m[G0.N+1];
	double *tau = new double[G0.N+1];
	double *G0_p = new double[G0.N+1];
	double *G0_m = new double[G0.N+1];

	fft_fermion_radix2_omega2tau(G0_p, G_omega, beta, G0.N, a);
	G0_p[G0.N] = -a - G0_p[0];

	for(int i=0; i<=G0.N; i++){
		tau[i] = (double)i * beta / (double)G0.N;
		G0_m[G0.N-i] = -G0_p[i];
	}

	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);

	delete [] tau;
	delete [] G0_p;
	delete [] G0_m;
}
/*
void G0_init_fft(struct green_func_0 &G0, complex<double> *G_omega, double beta, double V_sqr)
{
	double tau[G0.N+1], G0_p[G0.N+1], G0_m[G0.N+1];

	fft_fermion_radix2_omega2tau(G0_p, G_omega, beta, G0.N);
	G0_p[G0.N] = -1.0 - G0_p[0];

	for(int i=0; i<=G0.N; i++)  G0_p[i] *= V_sqr;

	for(int i=0; i<=G0.N; i++){
		tau[i] = (double)i * beta / (double)G0.N;
		G0_m[G0.N-i] = -G0_p[i];
	}

	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);
}
*/

void G0_init_fft(struct green_func_0 &G0, complex<double> *G_omega, double beta)
{
	G0_init_fft(G0, G_omega, beta, 1.0);
}
/*
void G0_init_fft(struct green_func_0 &G0, complex<double> *G_omega, double beta)
{
	double tau[G0.N+1], G0_p[G0.N+1], G0_m[G0.N+1];

	fft_fermion_radix2_omega2tau(G0_p, G_omega, beta, G0.N);
	G0_p[G0.N] = -1.0 - G0_p[0];

	for(int i=0; i<=G0.N; i++){
		tau[i] = (double)i * beta / (double)G0.N;
		G0_m[G0.N-i] = -G0_p[i];
	}

	gsl_spline_init(G0.spline_p, tau, G0_p, G0.N+1);
	gsl_spline_init(G0.spline_m, tau, G0_m, G0.N+1);
}
*/

/*
double G0_calc_interp(struct green_func_0 &G0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	#if EQUAL_TIME_G==0
	if(dif>=0)  return( gsl_spline_eval(G0.spline_p, dif, G0.acc_p) );
	#else
	if(dif>0)  return( gsl_spline_eval(G0.spline_p, dif, G0.acc_p) );
	#endif
	else  return( gsl_spline_eval(G0.spline_m, -dif, G0.acc_m) );
}
*/

// return G0(tau2-tau1)
// if tau2-tau1==0, return G0(+0)
double G0_calc_interp(struct green_func_0 &G0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	if(dif>=0)  return( gsl_spline_eval(G0.spline_p, dif, G0.acc_p) );
	else  return( gsl_spline_eval(G0.spline_m, -dif, G0.acc_m) );
}
// if tau2-tau1==0, return G0(-0)
double G0_calc_interp_0m(struct green_func_0 &G0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	if(dif>0)  return( gsl_spline_eval(G0.spline_p, dif, G0.acc_p) );
	else  return( gsl_spline_eval(G0.spline_m, -dif, G0.acc_m) );
}



void D0_alloc(struct boson_func_0 &D0, int N)
{
	D0.N = N;
	D0.acc_p = gsl_interp_accel_alloc ();
	D0.acc_m = gsl_interp_accel_alloc ();
	D0.spline_p = gsl_spline_alloc (INTERP, D0.N+1);
	D0.spline_m = gsl_spline_alloc (INTERP, D0.N+1);
}

void D0_free(struct boson_func_0 &D0)
{
	gsl_spline_free (D0.spline_p);
	gsl_spline_free (D0.spline_m);
	gsl_interp_accel_free (D0.acc_p);
	gsl_interp_accel_free (D0.acc_m);
}

// a = iw * D(iw), w->inf
void D0_init_fft(struct boson_func_0 &D0, complex<double> *D_omega, double beta, double a)
{
	double tau[D0.N+1], D0_p[D0.N+1], D0_m[D0.N+1];

	fft_boson_radix2_omega2tau(D0_p, D_omega, beta, D0.N, a);
	D0_p[D0.N] = D0_p[0] + a;

	for(int i=0; i<=D0.N; i++){
		tau[i] = (double)i * beta / (double)D0.N;
		D0_m[D0.N-i] = D0_p[i];
	}

	gsl_spline_init(D0.spline_p, tau, D0_p, D0.N+1);
	gsl_spline_init(D0.spline_m, tau, D0_m, D0.N+1);
}
void D0_init_fft(struct boson_func_0 &D0, complex<double> *D_omega, double beta)
{
	D0_init_fft(D0, D_omega, beta, 0);
}

/*
double D0_calc_interp(struct boson_func_0 &D0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	#if EQUAL_TIME_G==0
	if(dif>=0)  return( gsl_spline_eval(D0.spline_p, dif, D0.acc_p) );
	#else
	if(dif>0)  return( gsl_spline_eval(D0.spline_p, dif, D0.acc_p) );
	#endif
	else  return( gsl_spline_eval(D0.spline_m, -dif, D0.acc_m) );

// 	return( gsl_spline_eval(D0.spline, fabs(dif), D0.acc) );
}
*/

// return D0(tau2-tau1)
// if tau2-tau1==0, return D0(+0)
double D0_calc_interp(struct boson_func_0 &D0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	if(dif>=0)  return( gsl_spline_eval(D0.spline_p, dif, D0.acc_p) );
	else  return( gsl_spline_eval(D0.spline_m, -dif, D0.acc_m) );
}
// if tau2-tau1==0, return D0(-0)
double D0_calc_interp_0m(struct boson_func_0 &D0, double tau2, double tau1)
{
	double dif = tau2 - tau1;

	if(dif>0)  return( gsl_spline_eval(D0.spline_p, dif, D0.acc_p) );
	else  return( gsl_spline_eval(D0.spline_m, -dif, D0.acc_m) );
}

void K0_init_sum(struct boson_func_0 &K0, complex<double> *D_omega, double beta)
{
	double tau[K0.N+1];
	for(int i=0; i<=K0.N; i++){
		tau[i] = (double)i * beta / (double)K0.N;
	}

	double K_tau[K0.N+1];
	K_tau[0] = K_tau[K0.N] = 0;
	for(int i=1; i<=K0.N/2; i++){
		double beta_tau = beta - tau[i];
// 		double t1 = beta/2. - tau[i];
// 		double t2 = beta/2.;

		K_tau[i] = -0.5 * real(D_omega[0]) * tau[i] * beta_tau;

// 		double s = -1;
// 		for(int j=1; j<=K0.N/2; j++){
// 			double omega_b = double(2*j) * M_PI / beta;
// 			K_tau[i] += -s * 2.* real(D_omega[j]) / (omega_b*omega_b) * ( cos(t1 * omega_b) - cos(t2 * omega_b) );
// 			s = -s;
// 		}

// 		for(int j=1; j<=K0.N/2; j++){
// 			double omega_b = double(2*j) * M_PI / beta;
// 			K_tau[i] += -real(D_omega[j]) / (omega_b*omega_b) * ( cos( tau * omega_b) + cos( beta_tau * omega_b) -2. );
// 		}

		for(int j=1; j<=K0.N/2; j++){
			double omega_b = double(2*j) * M_PI / beta;
			K_tau[i] += 2.* real(D_omega[j]) / (omega_b*omega_b) * ( 1.- cos(tau[i] * omega_b) );  // WARNING: ACCURACY
		}

		K_tau[i] /= beta;

		K_tau[K0.N-i] = K_tau[i];
	}
// 	FILE *fp;
// 	char filename[128];
// 	sprintf(filename, DATA_DIR "%02d-K0_tau2.dat", 0);
// 	fp=fopen(filename, "w");
// 	for(int i=0; i<=K0.N; i++){
// 		fprintf(fp, "%.6e %.6e", tau[i], K_tau[i]);
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf("\n '%s'\n", filename);

	// symmetric
	gsl_spline_init(K0.spline_p, tau, K_tau, K0.N+1);
	gsl_spline_init(K0.spline_m, tau, K_tau, K0.N+1);
}

void K1_init_sum(struct boson_func_0 &K1, complex<double> *D_omega, double beta)
{
	double tau[K1.N+1];
	for(int i=0; i<=K1.N; i++){
		tau[i] = (double)i * beta / (double)K1.N;
	}

	double K_tau[K1.N+1];
	for(int i=0; i<=K1.N/2; i++){
		K_tau[i] = -real(D_omega[0]) * (0.5 - tau[i] / beta);

		for(int j=1; j<=K1.N/2; j++){
			double omega_b = double(2*j) * M_PI / beta;
			K_tau[i] += 2.* real(D_omega[j]) * sin(tau[i] * omega_b) / (omega_b * beta);
		}

		K_tau[K1.N-i] = -K_tau[i];
	}

// 	FILE *fp;
// 	char filename[128];
// 	sprintf(filename, DATA_DIR "%02d-K1_tau2.dat", 0);
// 	fp=fopen(filename, "w");
// 	for(int i=0; i<=K1.N; i++){
// 		fprintf(fp, "%.6e %.6e", tau[i], K_tau[i]);
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf("\n '%s'\n", filename);

	double K_tau_m[K1.N+1];
	double K_tau0 = K_tau[0];
	for(int i=0; i<=K1.N; i++){
		K_tau[i] -= K_tau0;  // boundary condition : K_tau[0] = 0
		K_tau_m[i] = -K_tau[i];  // antisymmetric
	}
// 	for(int i=0; i<=K1.N; i++){
// 		K_tau_m[i] = K_tau[K1.N-i];
// 	}
	gsl_spline_init(K1.spline_p, tau, K_tau, K1.N+1);
	gsl_spline_init(K1.spline_m, tau, K_tau_m, K1.N+1);

// 	FILE *fp;
// 	char filename[128];
// 	sprintf(filename, DATA_DIR "%02d-K1_tau2.dat", 0);
// 	fp=fopen(filename, "w");
// 	for(int i=0; i<=K1.N; i++){
// 		fprintf(fp, "%.6e %.6e\n", tau[i]-beta, D0_calc_interp(K1, tau[i]-beta, 0));
// 	}
// 	for(int i=0; i<=K1.N; i++){
// 		fprintf(fp, "%.6e %.6e\n", tau[i], D0_calc_interp(K1, tau[i], 0));
// 	}
// 	fclose(fp);
// 	printf("\n '%s'\n", filename);
}
