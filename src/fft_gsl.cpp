/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "fft.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>

using namespace std;

static const complex<double> IMAG(0, 1.0);

//
// G_tau   : real number, N data
// G_omega : complex number, N/2 data
//

static void fft_fermion_radix2_tau2omega_sub(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	double *G_temp = new double[4*N];
	
	//
	// extend range [0:beta) to [0:2*beta) with use of the anti-periodicity
	//
	for(int i=0; i<N; i++){
		G_temp[2*i] = G_tau[i];  // real part  [0:N)
		G_temp[2*i+1] = 0;  // imaginary part  [0:N)
		
		G_temp[2*(i+N)] = -G_tau[i];  // real part  [N:2N)
		G_temp[2*(i+N)+1] = 0;  // imaginary part  [N:2N)
	}
	
	gsl_fft_complex_radix2_backward(G_temp, 1, 2*N);
	
	//
	// pick up fermion Matsubara frequencies
	//
	double c = 0.5 * beta / (double)N;
	for(int i=0; i<N/2; i++){
		complex<double> temp(G_temp[2*(2*i+1)], G_temp[2*(2*i+1)+1]);
		G_omega[i] = c * temp;
	}
	
	delete [] G_temp;
}

//
// perform FFT after subtracting discontinuous function
// 
// a = G(-0) - G(+0)
//
void fft_fermion_radix2_tau2omega(double *G_tau, complex<double> *G_omega, double beta, int N, double a)
{
	double *G_tau_dif = new double[N];
	
// 	double e0 = log( -G_tau[0] / (a + G_tau[0]) ) / beta;
	
	double delta_tau = beta / (double)N;
	for(int i=0; i<N; i++){
		double tau = delta_tau * (double)i;
// 		G_tau_dif[i] = G_tau[i] + a * exp(-e0*tau) / ( exp(-e0*beta) + 1.0 );
		G_tau_dif[i] = G_tau[i] + a * 0.5;
	}
	
	fft_fermion_radix2_tau2omega_sub(G_tau_dif, G_omega, beta, N);
	
	for(int i=0; i<N/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
// 		G_omega[i] += a / (IMAG*omega - e0);
		G_omega[i] += a / (IMAG*omega);
	}
	
	delete [] G_tau_dif;
}

void fft_fermion_radix2_tau2omega(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	fft_fermion_radix2_tau2omega(G_tau, G_omega, beta, N, 1.0);
}


static void fft_fermion_radix2_omega2tau_sub(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	double *G_temp = new double[4*N];
	
	for(int i=0; i<N/2; i++){
		G_temp[4*i] = 0;  // real part  (boson)
		G_temp[4*i+1] = 0;  // imaginary part  (boson)
		
		G_temp[4*i+2] = real(G_omega[i]);  // real part (fermion)
		G_temp[4*i+3] = imag(G_omega[i]);  // imaginary part (fermion)
		
		//
		// negative frequencies
		//
		G_temp[4*(i+N/2)] = 0;  // real part  (boson)
		G_temp[4*(i+N/2)+1] = 0;  // imaginary part  (boson)
		
		G_temp[4*(i+N/2)+2] = real(G_omega[N/2-1-i]);  // real part (fermion)
		G_temp[4*(i+N/2)+3] = -imag(G_omega[N/2-1-i]);  // imaginary part (fermion)
	}
	
	gsl_fft_complex_radix2_forward(G_temp, 1, 2*N);
	
	//
	// extract real data on [0:beta)
	//
	double t = 1.0 / beta;
	for(int i=0; i<N; i++){
		G_tau[i] = t * G_temp[2*i];
	}
	
	delete [] G_temp;
}

//
// perform FFT after subtracting discontinuous function
// 
// a = G(-0) - G(+0)
//
void fft_fermion_radix2_omega2tau(double *G_tau, complex<double> *G_omega, double beta, int N, double a)
{
	complex<double> *G_omega_dif = new complex<double>[N/2];
	
	for(int i=0; i<N/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		G_omega_dif[i] = G_omega[i] - a / (IMAG*omega);
	}
	
	fft_fermion_radix2_omega2tau_sub(G_tau, G_omega_dif, beta, N);
	
	for(int i=0; i<N; i++){
		G_tau[i] -= 0.5 * a;
	}
	
	delete [] G_omega_dif;
}

// a = 1
void fft_fermion_radix2_omega2tau(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	fft_fermion_radix2_omega2tau(G_tau, G_omega, beta, N, 1.0);
}


/*
void fft_fermion_radix2_omega2tau(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	double G_temp1[4*N], G_temp2[4*N];
	
	for(int i=0; i<N; i++){
		G_temp1[4*i] = 0;  // real part  (boson)
		G_temp1[4*i+1] = 0;  // imaginary part  (boson)
		
		G_temp1[4*i+2] = real(G_omega[i]);  // real part (fermion)
		G_temp1[4*i+3] = imag(G_omega[i]);  // imaginary part (fermion)
	}
	for(int i=0; i<N; i++){
		G_temp2[4*i] = 0;  // real part  (boson)
		G_temp2[4*i+1] = 0;  // imaginary part  (boson)
		
		G_temp2[4*i+2] = real(G_omega[i]);  // real part (fermion)
		G_temp2[4*i+3] = -imag(G_omega[i]);  // imaginary part (fermion)
	}
	
	gsl_fft_complex_radix2_forward(G_temp1, 1, 2*N);
	gsl_fft_complex_radix2_backward(G_temp2, 1, 2*N);
	
	//
	// extract real data on [0:beta)
	//
	double t = 1.0 / beta;
	for(int i=0; i<N; i++){
		G_tau[i] = t * (G_temp1[2*i] + G_temp2[2*i]);
	}
}
*/


//
// G_tau   : real number, N data
// G_omega : complex number, N/2+1 data
//

void fft_boson_radix2_tau2omega(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	double *G_temp = new double[2*N];
	
	for(int i=0; i<N; i++){
		G_temp[2*i] = G_tau[i];  // real part  [0:N)
		G_temp[2*i+1] = 0;  // imaginary part  [0:N)
	}
	
	gsl_fft_complex_radix2_backward(G_temp, 1, N);
	
	double c = beta / (double)N;
	for(int i=0; i<=N/2; i++){
		complex<double> temp(G_temp[2*i], G_temp[2*i+1]);
		G_omega[i] = c * temp;
	}
	
	delete [] G_temp;
}

// a = G(-0) - G(+0)
void fft_boson_radix2_tau2omega(double *G_tau, complex<double> *G_omega, double beta, int N, double a)
{
// 	double G_tau_dif[N];
	double *G_tau_dif = new double[N];
	
	double delta_tau = beta / (double)N;
	for(int i=0; i<N; i++){
		double tau = delta_tau * (double)i;
// 		G_tau_dif[i] = G_tau[i] - a * exp(-tau / beta) / ( exp(-1) - 1.0 );
		G_tau_dif[i] = G_tau[i] - a * tau / beta;
	}
	
	fft_boson_radix2_tau2omega(G_tau_dif, G_omega, beta, N);
	
// 	for(int i=0; i<=N/2; i++){
// 		double omega = (double)(2*i) * M_PI / beta;
// 		G_omega[i] += a / (IMAG*omega - 1./ beta);
// 	}
	
	G_omega[0] += a * beta / 2.;
	for(int i=1; i<=N/2; i++){
		double omega = (double)(2*i) * M_PI / beta;
		G_omega[i] += a / (IMAG*omega);
	}
	
// 	FILE *fp;
// 	char filename[128];
// 	
// 	sprintf(filename, "data/" "%02d-test.dat", 0);
// 	fp=fopen(filename, "w");
// 	for(int i=0; i<N; i++){
// 		fprintf(fp, "%d", i);
// 		fprintf(fp, " %.5e", G_tau[i]);
// 		fprintf(fp, " %.5e", G_tau_dif[i]);
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	delete [] G_tau_dif;
}

void fft_boson_radix2_omega2tau(double *G_tau, complex<double> *G_omega, double beta, int N)
{
	double *G_temp = new double[2*N];
	
// 	for(int i=0; i<N/2; i++){
// 		G_temp[2*i] = real(G_omega[i]);  // real part  (boson)
// 		G_temp[2*i+1] = imag(G_omega[i]);  // imaginary part  (boson)
// 		
// 		//
// 		// negative frequencies
// 		//
// 		G_temp[2*(i+N/2)] = real(G_omega[N/2-1-i]);  // real part  (boson)
// 		G_temp[2*(i+N/2)+1] = -imag(G_omega[N/2-1-i]);  // imaginary part  (boson)
// 	}
	
	
	for(int i=0; i<=N/2; i++){
		G_temp[2*i] = real(G_omega[i]);  // real part  (boson)
		G_temp[2*i+1] = imag(G_omega[i]);  // imaginary part  (boson)
	}
	//
	// negative frequencies
	//
	for(int i=1; i<N/2; i++){
		G_temp[2*(N-i)] = real(G_omega[i]);  // real part  (boson)
		G_temp[2*(N-i)+1] = -imag(G_omega[i]);  // imaginary part  (boson)
	}
	
	gsl_fft_complex_radix2_forward(G_temp, 1, N);
	
	//
	// extract real data on [0:beta)
	//
	double t = 1.0 / beta;
	for(int i=0; i<N; i++){
		G_tau[i] = t * G_temp[2*i];
	}
	
	delete [] G_temp;
}

// a = G(-0) - G(+0)
void fft_boson_radix2_omega2tau(double *G_tau, complex<double> *G_omega, double beta, int N, double a)
{
	complex<double> *G_omega_dif = new complex<double>[N/2+1];
	
	G_omega_dif[0] = G_omega[0] - a * beta / 2.;
	
	for(int i=1; i<=N/2; i++){
		double omega = (double)(2*i) * M_PI / beta;
		G_omega_dif[i] = G_omega[i] - a / (IMAG*omega);
	}
	
	fft_boson_radix2_omega2tau(G_tau, G_omega_dif, beta, N);
	
	double delta_tau = beta / (double)N;
	for(int i=0; i<N; i++){
		double tau = delta_tau * (double)i;
		G_tau[i] += a * tau / beta;
	}
	
	delete [] G_omega_dif;
}
