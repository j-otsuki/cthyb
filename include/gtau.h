/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _GTAU_H
#define _GTAU_H

// #include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <complex>


class GTau{
public:
	// N+1
	void init_gtau(const std::vector<double>& g_tau, const std::vector<double>& tau);

	// N/2
	// a = iw * G(iw), w->inf
	void init_giw(const std::vector<std::complex<double> >& g_iw, double beta, double a);

	// return G(tau2, tau1)
	double calc_interp(double tau2, double tau1);

	~GTau();

private:
	gsl_interp_accel *acc_p;
	gsl_interp_accel *acc_m;
	gsl_spline *spline_p;  // G0(tau)  [0:beta]
	gsl_spline *spline_m;  // G0(-tau)  [0:beta]
	bool allocated = false;

	void alloc(std::size_t N);
	void free();
};


// interporation algorithm
#define INTERP gsl_interp_linear
#define INTERP_CHR "gsl_interp_linear"

// #define INTERP gsl_interp_cspline
// #define INTERP_CHR "gsl_itnerp_cspline"

// #define INTERP gsl_interp_akima
// #define INTERP_CHR "gsl_interp_akima"


#endif // _GTAU_H
