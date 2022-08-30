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
#include <cassert>
#include "fft.h"


using namespace std;

GTau::~GTau()
{
	free();
}

void GTau::alloc(size_t N)
{
	// G0.N = N;
	acc_p= gsl_interp_accel_alloc ();
	acc_m= gsl_interp_accel_alloc ();
	spline_p = gsl_spline_alloc (INTERP, N+1);  // Delta(tau)  [0:beta]
	spline_m = gsl_spline_alloc (INTERP, N+1);  // Delta(-tau)  [0:beta]
	allocated = true;
}

void GTau::free()
{
	if(allocated){
		gsl_spline_free (spline_p);
		gsl_spline_free (spline_m);
		gsl_interp_accel_free (acc_p);
		gsl_interp_accel_free (acc_m);
	}
	allocated = false;
}

void GTau::init_gtau(const std::vector<double>& g_tau, const std::vector<double>& tau)
{
	size_t N = g_tau.size() - 1;
	double beta = g_tau.back();
	assert (tau.size() == N + 1);

	free();
	alloc(N);

	vector<double> g_tau_m(N);
	for(int i=0; i<=N; i++){
		g_tau_m[N-i] = -g_tau[i];
	}

	gsl_spline_init(spline_p, tau.data(), g_tau.data(), N+1);
	gsl_spline_init(spline_m, tau.data(), g_tau_m.data(), N+1);
}

// a = iw * G(iw), w->inf
void GTau::init_giw(const std::vector<complex<double> >& g_iw, double beta, double a)
{
	size_t N = g_iw.size() * 2;

	// auto g_iw_noconst = const_cast<vector<complex<double> >&>(g_iw);  // remove const

	// FFT; G(iw) -> G(tau)
	// vector<double> g_tau(N+1);
	// fft_fermion_radix2_omega2tau(g_tau.data(), g_iw_noconst.data(), beta, N, a);
	// g_tau[N] = -a - g_tau[0];
	vector<double> g_tau;
	fft_fermion_iw2tau(g_tau, g_iw, beta, a);
	assert (g_tau.size() == N);
	g_tau.push_back(-a - g_tau[0]);

	vector<double> tau(N+1);
	for(int i=0; i<=N; i++){
		tau[i] = (double)i * beta / (double)N;
	}

	init_gtau(g_tau, tau);
}

// return G0(tau2-tau1)
// if tau2-tau1==0, return G0(+0)
double GTau::calc_interp(double tau2, double tau1)
{
	double dif = tau2 - tau1;

	if(dif>=0)  return( gsl_spline_eval(spline_p, dif, acc_p) );
	else  return( gsl_spline_eval(spline_m, -dif, acc_m) );
}
