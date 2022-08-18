/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "green_func_0.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "erfc.h"
#include "dmft.h"  // for Hilbert transform


// tau >= 0
double G0_calc_integ(double tau, double beta, double D)
{
// 	printf("\nInitializing G0 with gsl_integration_qag(qagi)\n");
// 	
	double result, error;
	double G0_func(double, void *);
	
	gsl_function F;
	F.function = &G0_func;
	double param[2] = {tau, beta};
	F.params = (void *)param;
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	
	gsl_integration_qag(&F, -D, D, 0, 1e-7, 1000, 6, w, &result, &error);
// 	gsl_integration_qagi(&F, 0, 1e-7, 1000, w, &result, &error);
	
	gsl_integration_workspace_free (w);
	
	return( result * 0.5 / D );
}

//
// integrand
//
double G0_func(double e, void *params)
{
// 	double tau = *(double *) params;
	
	double *param = (double *) params;
	double tau = param[0];
	double beta = param[1];
	
// 	if(tau >= 0){
		if(e >= 0){
			return( -exp(-e*tau) / ( exp(-e*beta)+1.0) );
		}
		else{
			return( -exp(e*(beta-tau)) / ( exp(e*beta)+1.0) );
		}
// 	}
// 	if(tau_i < 0){
// 		if(z >= 0){
// 			return( -dos_0 * exp(z*(prm.beta+tau_i)) / ( exp(-z*prm.beta)+1.0) );
// 		}
// 		else{
// 			return( -dos_0 * exp(-z*tau_i) / ( exp(z*prm.beta)+1.0) );
// 		}
// 	}
	
// 	double dos = sqrt(pow(prm.D, 2) - pow(e, 2));
// 	if(e >= 0){
// 		return( -exp(-e*tau) / ( exp(-e*prm.beta)+1.0) * dos );
// 	}
// 	else{
// 		return( -exp(e*(prm.beta-tau)) / ( exp(e*prm.beta)+1.0) * dos );
// 	}
}

void G0_calc_fft(double *G0_tau, double beta, double D, double V_sqr)
{
	complex<double> G0_omega[N_TAU/2];
	
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		
		G0_omega[i] = -IMAG * atan(D / omega) / D;  // = g_c
	}
	
	fft_fermion_radix2_omega2tau(G0_tau, G0_omega, beta, N_TAU);
	
	for(int i=0; i<N_TAU; i++){
		G0_tau[i] *= V_sqr;
	}
}



inline static void G0_omega_calc0(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	//
	// rectangular model
	//
// 	for(int i=0; i<N_TAU/2; i++){
// 		double omega = (double)(2*i+1) * M_PI / beta;
// 		G0_omega[i] = -IMAG / D * atan(D/omega);
// 	}
	
	//
	// rectangular model
	//
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		double x1 = D - chem_pot;
		double x2 = - D - chem_pot;
		double y = omega_n;
		
		double re_gc = - 0.25 * log( (x1*x1 + y*y) / (x2*x2 + y*y) ) / D;
		double im_gc = - 0.5 * ( atan(x1/y) - atan(x2/y) ) / D;
		
		G0_omega[i] = re_gc + IMAG * im_gc;
	}
}
inline static void G0_omega_calc1(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	//
	// semicircle
	//
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
// 		double x = omega / D;
		complex<double> z = (omega - IMAG * chem_pot) / D;
		G0_omega[i] = - 2.0 * IMAG / D * (sqrt(1.0 + z*z) - z);
	}
}
inline static void G0_omega_calc2(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	//
	// Gaussian
	//
	for(int i=0; i<N_TAU/2; i++){
// 		complex<double> i_omega(0, (double)(2*i+1) * M_PI / beta);
		complex<double> z = IMAG * (double)(2*i+1) * M_PI / beta + chem_pot;
		double fac = sqrt(2) / D;
		G0_omega[i] = fac * erfc_w(fac * z);
	}
}

inline static void G0_omega_calc_default(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	for(int i=0; i<N_TAU/2; i++){
		complex<double> z = IMAG * (double)(2*i+1) * M_PI / beta + chem_pot;
		G0_omega[i] = hilbert_transform(z, D);
	}
}

inline static void G0_omega_calc_m1(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	//
	// single bath ( G0(iw) = 1/(iw+mu),  dos = delta(omega+mu) )
	//
	for(int i=0; i<N_TAU/2; i++){
		complex<double> z = IMAG * (double)(2*i+1) * M_PI / beta + chem_pot;
		G0_omega[i] = 1./ z;
	}
}


// return filling
double G0_site_diag(complex<double> *G0_omega, double beta, double chem_pot, double D)
{
	switch(DOS){
		case 0:
			G0_omega_calc0(G0_omega, beta, chem_pot, D);
			break;
		case 1:
			G0_omega_calc1(G0_omega, beta, chem_pot, D);
			break;
		case 2:
			G0_omega_calc2(G0_omega, beta, chem_pot, D);
			break;
		case -1:
			G0_omega_calc_m1(G0_omega, beta, chem_pot, D);
			break;
		default:
			G0_omega_calc_default(G0_omega, beta, chem_pot, D);
			break;
	}
	
	return( eval_G_tau0m(G0_omega, beta) );
}

struct params_G0_site_diag{
// 	complex<double> G0_omega[][N_TAU/2];
	complex<double> (*G0_omega)[N_TAU/2];
	int N;
	double beta;
	double D;
	double n_c;
	double *E_c;
};
static double func_G0_site_diag(double chem_pot, void *params)
{
	struct params_G0_site_diag *p = (params_G0_site_diag *)params;
	
	double temp_n_c = 0;
	for(int a=0; a<p->N; a++){
		temp_n_c += G0_site_diag(p->G0_omega[a], p->beta, chem_pot-p->E_c[a], p->D);
	}
	temp_n_c /= (double)p->N;
	
	return( p->n_c - temp_n_c );
}

void G0_omega_calc(complex<double> G0_omega[][N_TAU/2], int N, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c)
{
	if( !flag_ensemble ){  // chemical potential given
		
		n_c = 0;
		for(int a=0; a<N; a++){
			n_c += G0_site_diag(G0_omega[a], beta, chem_pot-E_c[a], D);
		}
		n_c /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // occupation given
		struct params_G0_site_diag params = {G0_omega, N, beta, D, n_c, E_c};
		
		gsl_function F;
		F.function = &func_G0_site_diag;
		F.params = &params;
		
		double x_lo = -10.* D, x_hi = 10.* D;
// 		double x_lo = -D, x_hi = D;
		gsl_root_fsolver *s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
		gsl_root_fsolver_set (s, &F, x_lo, x_hi);
		
// 		printf ("\nusing %s method\n", gsl_root_fsolver_name (s));
// 		printf ("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root", "err");
		
		int status, iter = 0;
		do{
			iter++;
			status = gsl_root_fsolver_iterate (s);
			chem_pot = gsl_root_fsolver_root (s);
			x_lo = gsl_root_fsolver_x_lower (s);
			x_hi = gsl_root_fsolver_x_upper (s);
			status = gsl_root_test_interval (x_lo, x_hi, 0, ACCURACY_FILLING);
			
// 			if (status == GSL_SUCCESS)  printf ("Converged:\n");
// 			printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, chem_pot, x_hi - x_lo);
		}
		while (status == GSL_CONTINUE && iter < FILLING_ITER_MAX);
		
// 		printf("%s", status);
		gsl_root_fsolver_free (s);
	}
}

/*
void G0_omega_calc(complex<double> G0_omega[][N_TAU/2], int N, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c)
{
	if( !flag_ensemble ){  // grand canonical ensemble
		
		n_c = 0;
		for(int a=0; a<N; a++){
			n_c += G0_site_diag(G0_omega[a], beta, chem_pot-E_c[a], D);
		}
		n_c /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = 0;
		for(int a=0; a<N; a++){
			temp_n_c += G0_site_diag(G0_omega[a], beta, chem_pot-E_c[a], D);
		}
		temp_n_c /= (double)N;
		
		double delta_n_c = temp_n_c - n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING ){
			
			printf(" temp_n_c=%.4e  chem_pot=%.4e\n", temp_n_c, chem_pot);
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = 0;
			for(int a=0; a<N; a++){
				temp_n_c += G0_site_diag(G0_omega[a], beta, chem_pot-E_c[a], D);
			}
			temp_n_c /= (double)N;
			
			delta_n_c = temp_n_c - n_c;
			
			iter++;
// 			if(iter>100) exit(0);
		}
		
// 		printf("\nChemical potential : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf(" Filling : %.5lf\n", temp_n_c);
	}
}
*/
/*
void G0_omega_calc(complex<double> *G0_omega, double beta, double D, int flag_dmft, double &n_c, double &chem_pot)
{
	if( flag_dmft%2 == 0 ){
		
		n_c = G0_site_diag(G0_omega, beta, chem_pot, D);
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = G0_site_diag(G0_omega, beta, chem_pot, D);
		double delta_n_c = temp_n_c - n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING ){
			
// 			printf(" temp_n_c=%.4e  chem_pot=%.4e\n", temp_n_c, chem_pot);
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = G0_site_diag(G0_omega, beta, chem_pot, D);
			delta_n_c = temp_n_c - n_c;
			
			iter++;
// 			if(iter>100) exit(0);
		}
		
// 		printf("\nChemical potential : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf(" Filling : %.5lf\n", temp_n_c);
	}
}
*/

void G0_omega_calc(complex<double> *G0_omega, double beta, double D)
{
	double n_c, chem_pot=0;
	
	G0_site_diag(G0_omega, beta, chem_pot, D);
}


void G0_omega_calc(complex<double> *G0_omega, double beta, double D, double V_sqr)
{
	G0_omega_calc(G0_omega, beta, D);
	
	for(int i=0; i<N_TAU/2; i++)  G0_omega[i] *= V_sqr;
}

//=========================================================================================================

// unperturbed G_f
// for DOS==0 && symmetric condition
void Gf_calc_fft(double *Gf_tau, double beta, double D, double epsilon_f, double U, double V_sqr)
{
	complex<double> Gf_omega[N_TAU/2];
	
	double ef = epsilon_f + U * 0.5;  // add the Hartree term
	
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		
		complex<double> Delta_omega = -IMAG * V_sqr * atan(D / omega) / D;
		
		Gf_omega[i] = 1.0 / (IMAG*omega - ef - Delta_omega);
	}
	
	fft_fermion_radix2_omega2tau(Gf_tau, Gf_omega, beta, N_TAU);
}

/*
double calc_Gf(double tau)
{
// 	printf("\nInitializing G0 with gsl_integration_qag(qagi)\n");
// 	
	double result, error;
	double func_Gf(double, void *);
	
	gsl_function F;
	F.function = &func_Gf;
	F.params = &tau;
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	
	gsl_integration_qag(&F, -prm.D, prm.D, 0, 1e-7, 1000, 6, w, &result, &error);
// 	gsl_integration_qagi(&F, 0, 1e-7, 1000, w, &result, &error);
	
	gsl_integration_workspace_free (w);
	
	return( result );
}
//
// integrand
//
double func_Gf(double e, void *params)
{
	double tau = *(double *) params;
	
// 	double *param = (double *) params;
// 	double tau = param[0];
// 	double spin = param[1];
	
	double ef = prm.epsilon_f + prm.U * 0.5;  // add the Hartree term
	double dos = prm.Delta0 / ( pow(e - ef - prm.Delta0 / M_PI * log(fabs( (e+prm.D)/(e-prm.D) )), 2)
	 + pow(prm.Delta0, 2) ) / M_PI;
	
// 	if(tau >= 0){
		if(e >= 0){
			return( -exp(-e*tau) / ( exp(-e*prm.beta)+1.0 ) * dos);
		}
		else{
			return( -exp(e*(prm.beta-tau)) / ( exp(e*prm.beta)+1.0 ) * dos);
		}
// 	}
// 	if(tau_i < 0){
// 		if(z >= 0){
// 			return( -dos_0 * exp(z*(prm.beta+tau_i)) / ( exp(-z*prm.beta)+1.0) );
// 		}
// 		else{
// 			return( -dos_0 * exp(-z*tau_i) / ( exp(z*prm.beta)+1.0) );
// 		}
// 	}
	
	return(0);
}
*/




static double func_D0_omega_calc1(double x, void *params)
{
	double *param = (double *)params;
	double gamma = param[0];
	double y = param[1];
	
	return( pow(x, gamma+1.) / ( y*y + x*x ) );
}

// DOS = delta(w - w0)
void D0_omega_calc0(complex<double> *D0_omega, double beta, double w0)
{
	for(int i=0; i<=N_TAU/2; i++){
		double omega_b = double(2*i) * M_PI / beta;
		D0_omega[i] = -2.* w0 / ( omega_b*omega_b + w0*w0 );
	}
}

// DOS \propto omega^{gamma}  (with cutoff)
void D0_omega_calc1(complex<double> *D0_omega, double beta, double cutoff, double gamma)
{
	D0_omega[0] = 1./ gamma;
	
	gsl_function F;
	F.function = &func_D0_omega_calc1;
	double params[2] = {gamma};
	F.params = params;
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	int key = 6;
	double result, error;
	
	for(int i=1; i<=N_TAU/2; i++){
		double omega_b = double(2*i) * M_PI / beta;
		params[1] = omega_b / cutoff;
		gsl_integration_qag (&F, 0, 1, 0, 1e-7, 1000, key, w, &result, &error);
		D0_omega[i] = result;
		
// 		if( gamma == 1.0 ){
// 			double exact = 1.0 - params[1] * atan(1./ params[1]);
// 			printf(" check D0_omega  %.4e\n", exact - result);
// 			if( fabs( exact - result ) > 1e-10 ){
// 				printf("*** gamma = 1\n");
// 				exit(0);
// 			}
// 		}
// 		if( gamma == 2.0 ){
// 			double y_sqr = params[1] * params[1];
// 			double exact = 0.5 - y_sqr * acoth( 1.+ 2.* y_sqr );
// 			printf(" check D0_omega  %.4e\n", exact - result);
// 			if( fabs( exact - result ) > 1e-10 ){
// 				printf("*** gamma = 1\n");
// 				exit(0);
// 			}
// 		}
	}
	gsl_integration_workspace_free (w);
	
	for(int i=0; i<=N_TAU/2; i++){
		D0_omega[i] *= -2.* (gamma+1.) / cutoff;
	}
	
}

