/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "dmft.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "erfc.h"
#include "pade.h"
#include <gsl/gsl_sf_erf.h>
#include "green_func_0.h"  // for G0_omega_calc
#include "common.h"

static const double GSL_EPSABS=1.0e-10;
static const double GSL_EPSREL=1.0e-5;
static const int GSL_KEY=6;

// Parameters for function 'tune_chem_pot'
static double MU_MIN = -10;
static double MU_MAX = 10;
static int PRINT_ITERATION = 0;
static double OCCUP_ACCURACY = 1e-6;
static int OCCUP_ITER_MAX = 100;


static double eval_chem_pot(double nc, double D);


// evaluate G(tau=-0), occupation number
double eval_G_tau0m(complex<double> *G_omega, double beta, double jump)
{
	// occupation number n(e)
	
	double s=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		s += real(G_omega[i]);
	}
	s *= 2.0 / beta;
	s += 0.5 * jump;
	
	return( s );
}

double eval_G_tau0m(complex<double> *G_omega, double beta)
{
	return( eval_G_tau0m(G_omega, beta, 1.0) );
}


//============================================================================


// Hilbert transformation (2D)
static complex<double> hilbert_transform_2d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			double ek = t[0] * 2.* ( cosk[i] + cosk[j] );
// 			double ek = t[0] * 2.* ( cosk[i] + cosk[j] )
// 			 + t[1] * 4. * cosk[i] * cosk[j];
			
			r += w[i] * w[j] / (z - ek);
		}
	}
	
	return( r / double(N_L*N_L) );
}

// Hilbert transformation (3D)
static complex<double> hilbert_transform_3d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] );
// 				double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] )
// 				 + t[1] * 8.* cosk[i] * cosk[j] * cosk[k];
				
				r += w[i] * w[j] * w[k] / (z - ek);
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L) );
}

// Hilbert transformation (4D)
static complex<double> hilbert_transform_4d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				for(int l=0; l<=N_L; l++){
					double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] );
// 					double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] )
// 					 + t[1] * 16.* cosk[i] * cosk[j] * cosk[k] * cosk[l];
					
					r += w[i] * w[j] * w[k] * w[l] / (z - ek);
				}
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L*N_L) );
}


static double dos(double x, double D)
{
	double rho;
	
	#if DOS==0
	rho = 0.5 / D;
	// -D < x < D
	
	#elif DOS==1
	rho = 2.0 / (M_PI*D*D) * sqrt(D*D - x*x);
	// -D < x < D
	
	#elif DOS==2
	rho = 1.0 / D * sqrt(2.0/M_PI) * exp(-2.0 * (x/D) * (x/D));
	// -inf < x < inf
	
	#endif
	
	return(rho);
}
static double func_hilbert_transform_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0];
	double D = real(a[1]);
	
	double f = real( 1./ ( z - x ) );
	
	return( f * dos(x, D) );
}
static double func_hilbert_transform_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0];
	double D = real(a[1]);
	
	double f = imag( 1./ ( z - x ) );
	
	return( f * dos(x, D) );
}

static complex<double> hilbert_transform_integ(complex<double> z, double D)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	complex<double> a[2] = {z, D};
	
	gsl_function F;
	F.params = a;
	
	double (*func[2])(double, void *) = {func_hilbert_transform_re, func_hilbert_transform_im};
	double result[2], error[2];
	
	for(int f=0; f<2; f++){
		F.function = func[f];
		
		#if DOS==0 || DOS==1
		gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
		
		#elif DOS==2
		gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
		
		#endif
	}
	
	gsl_integration_workspace_free(w);
	
	return( result[0] + result[1] * IMAG );
}

static complex<double> hilbert_transform_bethe(complex<double> z, double D)
{
	z /= D;
	
	complex<double> r = -IMAG * sqrt( 1.0 - z*z );
	if(imag(z) < 0)  r = -r;
	
	r += z;
	r *= 2.0 / D;
	
	return r ;
}

// Hilbert transform
// output:
//     (1/N) sum_k (z - ek)^{-1} = int dos(e) / (z - e)
// input: z
//        D (half bandwidth or n.n. hopping)
// 
complex<double> hilbert_transform(complex<double> z, double D)
{
	double t[2] = {D, 0};
	
	#if DOS==0 || DOS==2
// 	#if DOS==0 || DOS==1 || DOS==2
		return hilbert_transform_integ(z, D);
	#elif DOS==1
		return hilbert_transform_bethe(z, D);
	#elif DOS==3
		return hilbert_transform_2d(z, t);
	#elif DOS==4
		return hilbert_transform_3d(z, t);
	#elif DOS==5
		return hilbert_transform_4d(z, t);
	#else
		printf("\n***error  in hilbert_transform(): no function\n");
		exit(0);
	#endif
}


//============================================================================
// Hilbert-square transformation
//   (1/N) sum_k (z - ek)^{-2} = int dos(e) / (z - e)^2


// 2D
static complex<double> hilbert_sqr_transform_2d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			double ek = t[0] * 2.* ( cosk[i] + cosk[j] );
			
			r += w[i] * w[j] / ((z - ek) * (z - ek));
		}
	}
	
	return( r / double(N_L*N_L) );
}

// 3D
static complex<double> hilbert_sqr_transform_3d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] );
				
				r += w[i] * w[j] * w[k] / ((z - ek) * (z - ek));
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L) );
}

// 4D
static complex<double> hilbert_sqr_transform_4d(complex<double> z, double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				for(int l=0; l<=N_L; l++){
					double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] );
					
					r += w[i] * w[j] * w[k] * w[l] / ((z - ek) * (z - ek));
				}
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L*N_L) );
}

static double func_hilbert_sqr_transform_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0];
	double D = real(a[1]);
	
	double f = real( 1./ (( z - x ) * ( z - x )) );
	
	return( f * dos(x, D) );
}
static double func_hilbert_sqr_transform_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0];
	double D = real(a[1]);
	
	double f = imag( 1./ (( z - x ) * ( z - x )) );
	
	return( f * dos(x, D) );
}

static complex<double> hilbert_sqr_transform_integ(complex<double> z, double D)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	complex<double> a[2] = {z, D};
	
	gsl_function F;
	F.params = a;
	
	double (*func[2])(double, void *) = {func_hilbert_sqr_transform_re, func_hilbert_sqr_transform_im};
	double result[2], error[2];
	
	for(int f=0; f<2; f++){
		F.function = func[f];
		
		#if DOS==0 || DOS==1
		gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
		
		#elif DOS==2
		gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
		
		#endif
	}
	
	gsl_integration_workspace_free(w);
	
	return( result[0] + result[1] * IMAG );
}

static complex<double> hilbert_sqr_transform_bethe(complex<double> z, double D)
{
	z /= D;
	
	complex<double> r = -hilbert_transform(z, D) / sqrt( 1.0 - z*z );
	r *= IMAG / D;
	if(imag(z) < 0)  r = -r;
	
	return r;
}

static complex<double> hilbert_sqr_transform_hypercube(complex<double> z, double D)
{
	return 4.0 / D / D * (z * hilbert_transform(z, D) - 1.0);
}



// Hilbert-sqr transform
// output:
//     (1/N) sum_k (z - ek)^{-2} = int dos(e) / (z - e)^2
// input: z
//        D (half bandwidth or n.n. hopping)
// 
complex<double> hilbert_sqr_transform(complex<double> z, double D)
{
	double t[2] = {D, 0};
	
	#if DOS==0
// 	#if DOS==0 || DOS==1 || DOS==2
		return -hilbert_sqr_transform_integ(z, D);
	#elif DOS==1
		return -hilbert_sqr_transform_bethe(z, D);
	#elif DOS==2
		return -hilbert_sqr_transform_hypercube(z, D);
	#elif DOS==3
		return -hilbert_sqr_transform_2d(z, t);
	#elif DOS==4
		return -hilbert_sqr_transform_3d(z, t);
	#elif DOS==5
		return -hilbert_sqr_transform_4d(z, t);
	#else
		printf("\n***error  in hilbert_sqr_transform(): no function\n");
		exit(0);
	#endif
}

//============================================================================

// Tune the chemical potential
// Return delta_mu
//
// INPUT:
//   double (*func_occup)(double delta_mu, void *params)
//     returns ( occup_calc - occup_given )
//
//   void *params
//     A parameter set used in func_occup
//
static double tune_chem_pot(double (*func_occup)(double, void *), void *params)
{
	double delta_mu;
	{
		gsl_function F;
		F.function = func_occup;
		F.params = params;
		
		double x_lo = MU_MIN, x_hi = MU_MAX;
		gsl_root_fsolver *s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
		gsl_root_fsolver_set (s, &F, x_lo, x_hi);
		
		if( PRINT_ITERATION>=2 ){
			printf ("\nusing %s method\n", gsl_root_fsolver_name (s));
			printf ("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root", "err");
		}
		
		int status, iter = 0;
		do{
			iter++;
			status = gsl_root_fsolver_iterate (s);
			delta_mu = gsl_root_fsolver_root (s);
			x_lo = gsl_root_fsolver_x_lower (s);
			x_hi = gsl_root_fsolver_x_upper (s);
			status = gsl_root_test_interval (x_lo, x_hi, 0, OCCUP_ACCURACY);
			
// 			if (status == GSL_SUCCESS)  printf ("Converged:\n");
			if( PRINT_ITERATION>=2 ){
				printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, delta_mu, x_hi - x_lo);
			}
		}
		while (status == GSL_CONTINUE && iter < OCCUP_ITER_MAX);
// 		printf("%s\n", status);
		
		gsl_root_fsolver_free (s);
	}
	return delta_mu;
}

//============================================================================


static inline complex<double> Gc_site_diag0(complex<double> zeta, double D)
{
	//
	// rectangular model
	//
// 	double x1 = D + real(self_c);
// 	double x2 = - D + real(self_c);
// 	double y = omega_n - imag(self_c);
	
	
// 	double x1 = D - real(zeta);
// 	double x2 = - D - real(zeta);
// 	double y = imag(zeta);
// 	
// 	double re_gc = - 0.25 * log( (x1*x1 + y*y) / (x2*x2 + y*y) ) / D;
// 	double im_gc = - 0.5 * ( atan(x1/y) - atan(x2/y) ) / D;
// 	
// 	return(re_gc + IMAG * im_gc);
	
	
	double x1 = real(zeta) + D;
	double x2 = real(zeta) - D;
	double y = fabs(imag(zeta));
	
	double re_gc = 0.25 * log( (x1*x1 + y*y) / (x2*x2 + y*y) ) / D;
	double im_gc = - 0.5 * ( atan(x1/y) - atan(x2/y) ) / D;
	
	if(imag(zeta) < 0)  im_gc = -im_gc;
	
	return(re_gc + IMAG * im_gc);
}

static inline complex<double> Gc_site_diag1(complex<double> zeta, double D)
{
	//
	// semicircle
	//
/*	complex<double> z1 = (IMAG * omega_n - self_c) / D;
	
// 	return( 2.0 / D * (z1 - sqrt(z1*z1 - 1.0)) );
	complex<double> r1 = 2.0 / D * (z1 - sqrt(z1*z1 - 1.0));
	complex<double> temp1 = sqrt(z1*z1 - 1.0);
	
	complex<double> z = (omega_n + IMAG * self_c) / D;
	
// 	return( -2.0 * IMAG / D * (sqrt(z*z + 1.0) - z) );
	complex<double> r2 = -2.0 * IMAG / D * (sqrt(z*z + 1.0) - z);
	complex<double> temp2 = IMAG * sqrt(z*z + 1.0);
	
	printf("(%.1e, %.1e) (%.1e, %.1e) ", real(r1), imag(r1), real(r2), imag(r2));
	printf("(%.1e, %.1e) (%.1e, %.1e)\n", real(temp1), imag(temp1), real(temp2), imag(temp2));
	return(r2);
*/	
	
// 	complex<double> z = (omega_n + IMAG * self_c) / D;
	
// 	complex<double> z = -IMAG * zeta / D;
// 	
// 	return( -2.0 * IMAG / D * (sqrt(z*z + 1.0) - z) );
	
	
	complex<double> z = zeta / D;
	
	complex<double> gc = -IMAG * sqrt( 1.0 - z*z );
	if(imag(zeta) < 0)  gc = -gc;
	gc += z;
	gc *= 2.0 / D;
	
	return(gc);
}

static inline complex<double> Gc_site_diag2(complex<double> zeta, double D)
{
	//
	// Gaussian
	//
// 	complex<double> z = (IMAG * omega_n - self_c);
	
	double fac = sqrt(2) / D;
	if(imag(zeta) < 0)  fac = -fac;
	
	return( fac * erfc_w(fac * zeta) );
}

//
// return band filling
//
static double Gc_site_diag(complex<double> *Gc_perio, complex<double> *self_c, double beta, double chem_pot, double D)
{
	
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		complex<double> zeta = IMAG * omega_n - self_c[i] + chem_pot;
		
		#if DOS==0
// 		Gc_perio[i] = Gc_site_diag0(omega_n, self_c[i] - chem_pot, D);
		Gc_perio[i] = Gc_site_diag0(zeta, D);
		
		#elif DOS==1
// 		Gc_perio[i] = Gc_site_diag1(omega_n, self_c[i] - chem_pot, D);
		Gc_perio[i] = Gc_site_diag1(zeta, D);
		
		#elif DOS==2
// 		Gc_perio[i] = Gc_site_diag2(omega_n, self_c[i] - chem_pot, D);
		Gc_perio[i] = Gc_site_diag2(zeta, D);
		
		#endif
	}
	
	// evaluate filling
	double Gtau_temp[N_TAU];
	fft_fermion_radix2_omega2tau(Gtau_temp, Gc_perio, beta, N_TAU);
	
	return(1.0 + Gtau_temp[0]);
}

static double Gc_site_diag_2sub(complex<double> *Gc_perio[2], complex<double> *self_c[2], double beta, double chem_pot, double D, double *n_c[2])
{
	
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		complex<double> zeta_A = IMAG * omega_n - self_c[0][i] + chem_pot;
		complex<double> zeta_B = IMAG * omega_n - self_c[1][i] + chem_pot;
		complex<double> zeta = sqrt(zeta_A * zeta_B);
// 		complex<double> zeta = -sqrt(zeta_A * zeta_B);
		
		complex<double> temp_Gc_perio;
		
		#if DOS==0
		temp_Gc_perio = Gc_site_diag0(zeta, D);
// 		temp_Gc_perio = (Gc_site_diag0(zeta, D) - Gc_site_diag0(-zeta, D)) / 2.;
		
		#elif DOS==1
		temp_Gc_perio = Gc_site_diag1(zeta, D);
// 		temp_Gc_perio = (Gc_site_diag1(zeta, D) - Gc_site_diag1(-zeta, D)) / 2.;
		
		#elif DOS==2
		temp_Gc_perio = Gc_site_diag2(zeta, D);
// 		temp_Gc_perio = (Gc_site_diag2(zeta, D) - Gc_site_diag2(-zeta, D)) / 2.;
		
		#endif
		
		temp_Gc_perio /= zeta;
		Gc_perio[0][i] = temp_Gc_perio * zeta_B;
		Gc_perio[1][i] = temp_Gc_perio * zeta_A;
	}
	
	// evaluate filling
	double Gtau_temp[N_TAU];
	for(int l=0; l<2; l++){
		fft_fermion_radix2_omega2tau(Gtau_temp, Gc_perio[l], beta, N_TAU);
		*n_c[l] = 1.0 + Gtau_temp[0];
	}
	
	return( (*n_c[0] + *n_c[1]) / 2. );
}


//
// for the Kondo & Coqblin-Schrieffer model
//  t = V^2 G_f
//
/*
void dmft(struct dmft_green_func &G, complex<double> *tmat, double beta, double D, int flag_dmft, double &n_c, double &chem_pot)
{
	for(int i=0; i<N_TAU/2; i++){
		G.self_c[i] = tmat[i] / ( 1.0 + G.Gc_cav_omega[i] * tmat[i] );
	}
	
	
	if( flag_dmft%2 == 0 ){
		
		double Gc_site_diag(complex<double> *, complex<double> *, double, double, double);
		
		n_c = Gc_site_diag(G.Gc_perio, G.self_c, beta, chem_pot, D);
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = Gc_site_diag(G.Gc_perio, G.self_c, beta, chem_pot, D);
		double delta_n_c = temp_n_c - n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING ){
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = Gc_site_diag(G.Gc_perio, G.self_c, beta, chem_pot, D);
			delta_n_c = temp_n_c - n_c;
			
			iter++;
		}
		
// 		printf("\nChemical potential : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf(" Filling : %.5lf\n", temp_n_c);
	}
	
	
	
	for(int i=0; i<N_TAU/2; i++){
		
		G.tmat_perio[i] = G.self_c[i] + pow(G.self_c[i], 2) * G.Gc_perio[i];
		
		
// 		#if DOS==1 // bethe lattice
// 		G.Gc_cav_omega[i] = 1./ ( IMAG * omega_n - 0.25 * D * D * G.Gc_perio[i] );
// 		
// 		#else
		G.Gc_cav_omega[i] = G.Gc_perio[i] / ( 1.0 + G.Gc_perio[i] * G.self_c[i] );
		
// 		#endif
		
		
// 		complex<double> temp_cav = G.Gc_perio[i] / ( 1.0 + G.Gc_perio[i] * G.self_c[i] );
// 		double r = 0.05;
// 		G.Gc_cav_omega[i] = (1.0-r) * G.Gc_cav_omega[i] + r * temp_cav;
	}
	
}
*/

//
// for Kondo & CS model
//
void dmft_cs(struct dmft_green_func &G, complex<double> *tmat, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot)
{
	complex<double> *p_tmat[1];
	p_tmat[0] = tmat;
	
	double E_c[1] = {0};
	
	dmft_cs(&G, p_tmat, 1, beta, D, flag_ensemble, n_c, chem_pot, E_c);
}

void dmft_cs(struct dmft_green_func *G, complex<double> **tmat, int N, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c)
{
	for(int a=0; a<N; a++){
		for(int i=0; i<N_TAU/2; i++){
			G[a].self_c[i] = tmat[a][i] / ( 1.0 + G[a].Gc_cav_omega[i] * tmat[a][i] );
		}
	}
	
	
	if( !flag_ensemble ){  // grand canonical ensemble
		
		n_c = 0;
		for(int a=0; a<N; a++){
			G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
			n_c += G[a].n_c;
		}
		n_c /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double chem_pot_old = chem_pot;
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = 0;
		for(int a=0; a<N; a++){
			G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
			temp_n_c += G[a].n_c;
		}
		temp_n_c /= (double)N;
		
		double delta_n_c = temp_n_c - n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING && iter < FILLING_ITER_MAX){
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = 0;
			for(int a=0; a<N; a++){
				G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
				temp_n_c += G[a].n_c;
			}
			temp_n_c /= (double)N;
			
			delta_n_c = temp_n_c - n_c;
			
			iter++;
		}
		
// 		printf("\nChemical potential changed : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf("Filling : %.5lf\n", temp_n_c);
		
/*		
		double CHEM_UPDATE_RATIO=1.0;
		chem_pot = CHEM_UPDATE_RATIO * chem_pot + (1.0-CHEM_UPDATE_RATIO) * chem_pot_old;
		
		temp_n_c = 0;
		for(int a=0; a<N; a++){
			G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
			temp_n_c += G[a].n_c;
		}
		temp_n_c /= (double)N;
		
		printf("\nChemical potential re-changed : %.5e  (%d iteration)\n", chem_pot, iter);
		printf("Filling : %.5lf\n", temp_n_c);
*/		
	}
	
	
	for(int a=0; a<N; a++){
		for(int i=0; i<N_TAU/2; i++){
			
// 			G[a].tmat_perio[i] = G[a].self_c[i] + pow(G[a].self_c[i], 2) * G[a].Gc_perio[i];
			G[a].tmat_perio[i] = G[a].self_c[i] * (1.0 + G[a].self_c[i] * G[a].Gc_perio[i]);
			
			
			#if DOS==1 // bethe lattice
			double omega_n = (double)(2*i+1) * M_PI / beta;
// 			G[a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n - 0.25 * D * D * G[a].Gc_perio[i] );
			G[a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n + chem_pot-E_c[a] - 0.25 * D * D * G[a].Gc_perio[i] );
			
			#else
			G[a].Gc_cav_omega[i] = G[a].Gc_perio[i] / ( 1.0 + G[a].Gc_perio[i] * G[a].self_c[i] );
			
			#endif
			
			
	// 		complex<double> temp_cav = G.Gc_perio[i] / ( 1.0 + G.Gc_perio[i] * G.self_c[i] );
	// 		double r = 0.05;
	// 		G.Gc_cav_omega[i] = (1.0-r) * G.Gc_cav_omega[i] + r * temp_cav;
		}
	}
	
}



void dmft_cs_2sub(struct dmft_green_func *G[2], complex<double> **tmat[2], int N, double beta, double D, int flag_ensemble, double &ave_n_c, double &chem_pot, double *E_c)
{
	for(int l=0; l<2; l++){
		for(int a=0; a<N; a++){
			for(int i=0; i<N_TAU/2; i++){
				G[l][a].self_c[i] = tmat[l][a][i] / ( 1.0 + G[l][a].Gc_cav_omega[i] * tmat[l][a][i] );
			}
		}
	}
	
	
	complex<double> *sub_Gc_perio[N][2], *sub_self_c[N][2];
	double *sub_nc[N][2];
	
	for(int a=0; a<N; a++){
		for(int l=0; l<2; l++){
			sub_Gc_perio[a][l] = G[l][a].Gc_perio;
			sub_self_c[a][l] = G[l][a].self_c;
			sub_nc[a][l] = &G[l][a].n_c;
		}
	}
	
	if( !flag_ensemble ){  // grand canonical ensemble
		
		ave_n_c = 0;
		for(int a=0; a<N; a++){
			ave_n_c += Gc_site_diag_2sub(sub_Gc_perio[a], sub_self_c[a], beta, chem_pot-E_c[a], D, sub_nc[a]);
		}
		ave_n_c /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double chem_pot_old = chem_pot;
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = 0;
		for(int a=0; a<N; a++){
			temp_n_c += Gc_site_diag_2sub(sub_Gc_perio[a], sub_self_c[a], beta, chem_pot-E_c[a], D, sub_nc[a]);
		}
		temp_n_c /= (double)N;
		
		double delta_n_c = temp_n_c - ave_n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING && iter < FILLING_ITER_MAX){
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = 0;
			for(int a=0; a<N; a++){
				temp_n_c += Gc_site_diag_2sub(sub_Gc_perio[a], sub_self_c[a], beta, chem_pot-E_c[a], D, sub_nc[a]);
			}
			temp_n_c /= (double)N;
			
			delta_n_c = temp_n_c - ave_n_c;
			
			iter++;
		}
		
// 		printf("\nChemical potential changed : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf("Filling : %.5lf\n", temp_n_c);
		
	}
	
	
	for(int l=0; l<2; l++){
		for(int a=0; a<N; a++){
			for(int i=0; i<N_TAU/2; i++){
				
				G[l][a].tmat_perio[i] = G[l][a].self_c[i] * (1.0 + G[l][a].self_c[i] * G[l][a].Gc_perio[i]);
				
				
				#if DOS==1 // bethe lattice
				double omega_n = (double)(2*i+1) * M_PI / beta;
	// 			G[a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n - 0.25 * D * D * G[a].Gc_perio[i] );
				G[l][a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n + chem_pot-E_c[a] - 0.25 * D * D * G[l][a].Gc_perio[i] );
				
				#else
				G[l][a].Gc_cav_omega[i] = G[l][a].Gc_perio[i] / ( 1.0 + G[l][a].Gc_perio[i] * G[l][a].self_c[i] );
				
				#endif
			}
		}
	}
	
}


//============================================================================
//
// for the Anderson lattice model (f-level fix)
// input: Gf
//
void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr)
{
	complex<double> tmat_omega[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  tmat_omega[i] = V_sqr * Gf[i];
	
	double temp_n_c, chem_pot=0;
	dmft_cs(G, tmat_omega, beta, D, 0, temp_n_c, chem_pot);
	
	for(int i=0; i<N_TAU/2; i++)  G.tmat_perio[i] /= V_sqr;
	
	
	// evaluate filling
	double Gf_tau_temp[N_TAU];
	fft_fermion_radix2_omega2tau(Gf_tau_temp, G.tmat_perio, beta, N_TAU);
	
	G.n_f = 1.0 + Gf_tau_temp[0];
	
// 	printf(" n_f=%.5lf  n_c=%.5lf\n", G.n_f, temp_n_c);
}

//
// n_f, tot_n: particle numbers per spin
//
void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot)
{
	dmft_ads0(G, Gf, beta, D, V_sqr, n_f, flag_ensemble, tot_n, chem_pot, 1.0);
}

void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double jump)
{
// 	G.Gf_perio = G.tmat_perio;
	
	complex<double> tmat_omega[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  tmat_omega[i] = V_sqr * Gf[i];
	
	double temp_n_c = tot_n - n_f;
	dmft_cs(G, tmat_omega, beta, D, flag_ensemble, temp_n_c, chem_pot);
	
	// void dmft_cs(struct dmft_green_func *G, complex<double> **tmat, int N,
	//              double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c)
	
	if( !flag_ensemble ){  // grand canonical ensemble
		tot_n = temp_n_c + n_f;
	}
	
	for(int i=0; i<N_TAU/2; i++)  G.tmat_perio[i] /= V_sqr;
	
	
	// evaluate filling
	double Gf_tau_temp[N_TAU];
	fft_fermion_radix2_omega2tau(Gf_tau_temp, G.tmat_perio, beta, N_TAU, jump);
	
// 	G.n_f = 1.0 + Gf_tau_temp[0];
	G.n_f = jump + Gf_tau_temp[0];
	
// 	printf("\n chem_pot=%.5e\n", chem_pot);
// 	printf(" tot_n=%.5lf  n_f=%.5lf  n_c=%.5lf\n", tot_n, G.n_f, G.n_c);
	
}


//
// n_f, tot_n: particle numbers (not per spin)
//
void dmft_ads0(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double *E_c, double *jump)
{
// 	for(int a=0; a<N; a++){
// 		G[a].Gf_perio = G[a].tmat_perio;
// 	}
	
	// void dmft_cs(struct dmft_green_func *G, complex<double> **tmat, int N,
	//              double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c)
	
	complex<double> *tmat_omega[N];
	for(int a=0; a<N; a++){
		tmat_omega[a] = new complex<double>[N_TAU/2];
		for(int i=0; i<N_TAU/2; i++)  tmat_omega[a][i] = V_sqr[a] * Gf[a][i];
	}
	
	double temp_n_c = (tot_n - n_f) / (double)N;
	dmft_cs(G, tmat_omega, N, beta, D, flag_ensemble, temp_n_c, chem_pot, E_c);
	
	if( !flag_ensemble ){  // grand canonical ensemble
		tot_n = n_f + temp_n_c * (double)N;
	}
	
	for(int a=0; a<N; a++){
		for(int i=0; i<N_TAU/2; i++)  G[a].tmat_perio[i] /= V_sqr[a];
	}
	
	// evaluate filling
	for(int a=0; a<N; a++){
		double Gf_tau_temp[N_TAU];
		fft_fermion_radix2_omega2tau(Gf_tau_temp, G[a].tmat_perio, beta, N_TAU, jump[a]);
		
	// 	G.n_f = 1.0 + Gf_tau_temp[0];
		G[a].n_f = jump[a] + Gf_tau_temp[0];
	}
	
	for(int a=0; a<N; a++){
		delete [] tmat_omega[a];
	}
	
// 	printf("\n chem_pot=%.5e\n", chem_pot);
// 	printf(" tot_n=%.5lf  n_f=%.5lf  n_c=%.5lf\n", tot_n, G.n_f, G.n_c);
}


//
// n_f, tot_n: particle numbers (not per spin)
//
void dmft_ads0_2sub(struct dmft_green_func *G[2], complex<double> **Gf[2], int N, double beta, double D, double *V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double *E_c, double *jump[2])
{
// 	for(int a=0; a<N; a++){
// 		G[0][a].Gf_perio = G[0][a].tmat_perio;
// 		G[1][a].Gf_perio = G[1][a].tmat_perio;
// 	}
	
	// void dmft_cs_2sub(struct dmft_green_func *G[2], complex<double> **tmat[2], int N,
	//              double beta, double D, int flag_ensemble, double &ave_n_c, double &chem_pot, double *E_c)
	
	complex<double> *tmat_A[N], *tmat_B[N];
	complex<double> **tmat[2] = {tmat_A, tmat_B};
	
	for(int l=0; l<2; l++){
		for(int a=0; a<N; a++){
			tmat[l][a] = new complex<double>[N_TAU/2];
			for(int i=0; i<N_TAU/2; i++)  tmat[l][a][i] = V_sqr[a] * Gf[l][a][i];
		}
	}
	
	double temp_n_c = (tot_n - n_f) / (double)N;
// 	dmft_cs(G, tmat, N, beta, D, flag_ensemble, temp_n_c, chem_pot, E_c);
	dmft_cs_2sub(G, tmat, N, beta, D, flag_ensemble, temp_n_c, chem_pot, E_c);
	
	if( !flag_ensemble ){  // grand canonical ensemble
		tot_n = n_f + temp_n_c * (double)N;
	}
	
	for(int l=0; l<2; l++){
		for(int a=0; a<N; a++){
			for(int i=0; i<N_TAU/2; i++)  G[l][a].tmat_perio[i] /= V_sqr[a];
		}
		
		// evaluate filling
		for(int a=0; a<N; a++){
			double Gf_tau_temp[N_TAU];
			fft_fermion_radix2_omega2tau(Gf_tau_temp, G[l][a].tmat_perio, beta, N_TAU, jump[l][a]);
			
		// 	G.n_f = 1.0 + Gf_tau_temp[0];
			G[l][a].n_f = jump[l][a] + Gf_tau_temp[0];
		}
		
		for(int a=0; a<N; a++)  delete [] tmat[l][a];
	}
	
	printf("\n chem_pot=%.5e\n", chem_pot);
	printf(" tot_n=%.5lf  n_f=%.5lf  n_c(per spin)=%.5lf\n", tot_n, n_f, temp_n_c);
}


/*
//
// n_f, tot_n: particle numbers (not per spin)
//
void dmft_ads0_2sublat(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double *E_c, double *jump)
{
	for(int a=0; a<N; a++){
		for(int i=0; i<N_TAU/2; i++){
			complex<double> tmat = V_sqr[a] * Gf[a][i];
			G[a].self_c[i] = tmat / ( 1.0 + tmat * G[a].Gc_cav_omega[i] );
		}
	}
	
	
	
	if( !flag_ensemble ){  // grand canonical ensemble
		
		double n_c = 0;
		for(int a=0; a<N; a++){
			G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
			n_c += G[a].n_c;
		}
// 		n_c /= (double)N;
		tot_n = n_f + n_c;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double n_c = tot_n - n_f;
		double chem_pot_old = chem_pot;
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n_c = 0;
		for(int a=0; a<N; a++){
			G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
			temp_n_c += G[a].n_c;
		}
		
		double delta_n_c = temp_n_c - n_c;
		
		int sign_delta = 1;
		if( delta_n_c < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n_c) > ACCURACY_FILLING && iter < FILLING_ITER_MAX){
			
			if( (double)sign_delta * delta_n_c < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n_c > 0 )  chem_pot -= d_chem_pot;
			else                 chem_pot += d_chem_pot;
			
			
			temp_n_c = 0;
			for(int a=0; a<N; a++){
				G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
				temp_n_c += G[a].n_c;
			}
			
			delta_n_c = temp_n_c - n_c;
			
			iter++;
		}
		
// 		printf("\nChemical potential changed : %.5e  (%d iteration)\n", chem_pot, iter);
// 		printf("Filling : %.5lf\n", temp_n_c);
	}
	
	
	
	for(int a=0; a<N; a++){
		for(int i=0; i<N_TAU/2; i++){
			
// 			G[a].tmat_perio[i] = G[a].self_c[i] + pow(G[a].self_c[i], 2) * G[a].Gc_perio[i];
			G[a].Gf_perio[i] = G[a].self_c[i] * (1.0 + G[a].self_c[i] * G[a].Gc_perio[i]) / V_sqr[a];
			
			
			#if DOS==1 // bethe lattice
			double omega_n = (double)(2*i+1) * M_PI / beta;
// 			G[a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n - 0.25 * D * D * G[a].Gc_perio[i] );
			G[a].Gc_cav_omega[i] = 1./ ( IMAG * omega_n + chem_pot-E_c[a] - 0.25 * D * D * G[a].Gc_perio[i] );
			
			#else
			G[a].Gc_cav_omega[i] = G[a].Gc_perio[i] / ( 1.0 + G[a].Gc_perio[i] * G[a].self_c[i] );
			
			#endif
		}
		
		// evaluate filling
		double Gf_tau_temp[N_TAU];
		fft_fermion_radix2_omega2tau(Gf_tau_temp, G[a].Gf_perio, beta, N_TAU, jump[a]);
		
// 		G.n_f = 1.0 + Gf_tau_temp[0];
		G[a].n_f = jump[a] + Gf_tau_temp[0];
	}
	
	
	
	
// 	printf("\n chem_pot=%.5e\n", chem_pot);
// 	printf(" tot_n=%.5lf  n_f=%.5lf  n_c=%.5lf\n", tot_n, G.n_f, G.n_c);
	
}
*/

//============================================================================
//
// for the Anderson lattice model
// dmft_ads1
//   ef -> ef - chem_pot
//   input: self_f
//

static double func_Gf_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = real( 1.0 / ( z2 - f_disp * x - V_sqr / (z1 - x - ec) ) );
	
	return( f * dos(x, D) );
}
static double func_Gf_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = imag( 1.0 / ( z2 - f_disp * x - V_sqr / (z1 - x - ec) ) );
	
	return( f * dos(x, D) );
}
static double func_Gc_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = real( 1.0 / ( z1 - x - ec - V_sqr / (z2 - f_disp * x) ) );
	
	return( f * dos(x, D) );
}
static double func_Gc_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = imag( 1.0 / ( z1 - x - ec - V_sqr / (z2 - f_disp * x) ) );
	
	return( f * dos(x, D) );
}


// return : Gf_perio[], Gc_perio[], n_f, n_c
static void G_site_diag_ads1(complex<double> *Gf_perio, complex<double> *Gc_perio, complex<double> *self_f, double &n_f, double &n_c, double beta, double D, double V_sqr, double ef, double f_disp, double ec, double chem_pot)
{
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	for(int i=0; i<N_TAU/2; i++){
// 		printf(" %d", i);
		
		complex<double> z1 = IMAG * (double)(2*i+1) * M_PI / beta;
		complex<double> z2 = z1 - ef + chem_pot - self_f[i];
// 		complex<double> z2 = z1 - ef - self_f[i];
		
		complex<double> a[6] = {z1, z2, V_sqr, f_disp, ec-chem_pot, D};
		
		gsl_function F;
		F.params = a;
		
		double (*func[4])(double, void *) = {func_Gf_re, func_Gf_im, func_Gc_re, func_Gc_im};
		double result[4], error[4];
		
// 		for(int j=0; j<201; j++){
// 			double x = -1.0 + (double)j * 0.01;
// 			printf("%d %.5lf %.5lf %.5lf %.5lf\n", j, func[0](x,a), func[1](x,a), func[2](x,a), func[3](x,a));
// 		}
		
		for(int f=0; f<4; f++){
// 			printf("%d", f);
			F.function = func[f];
			
			#if DOS==0 || DOS==1
// 			gsl_integration_qags(&F, -D, D, 0, 1e-7, 1000, w, &result[f], &error[f]);
			gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
			
			#elif DOS==2
			gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
			
			#endif
			
// 			printf("  %.6e %.6e", result[f], error[f]);
		}
// 		printf("\n");
		
		
		Gf_perio[i] = result[0] + result[1] * IMAG;
		Gc_perio[i] = result[2] + result[3] * IMAG;
	}
	
	gsl_integration_workspace_free(w);
	
	
	// evaluate filling
	n_f = eval_G_tau0m(Gf_perio, beta);
	n_c = eval_G_tau0m(Gc_perio, beta);
	
// 	printf(" n_f=%.5lf  n_c=%.5lf  mu=%.5lf\n", n_f, n_c, chem_pot);
	
}



static double Veff_2d(double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	double r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			double ek = t[0] * 2.* ( cosk[i] + cosk[j] )
			 + t[1] * 4. * cosk[i] * cosk[j];
			
			r += w[i] * w[j] * ek * ek;
		}
	}
	
	return( r / double(N_L*N_L) );
}

static double Veff_3d(double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	double r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] )
				 + t[1] * 8.* cosk[i] * cosk[j] * cosk[k];
				
				r += w[i] * w[j] * w[k] * ek * ek;
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L) );
}

static double Veff_4d(double t[2])
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	double r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				for(int l=0; l<=N_L; l++){
					double ek = t[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] );
					 + t[1] * 16.* cosk[i] * cosk[j] * cosk[k] * cosk[l];
					
					r += w[i] * w[j] * w[k] * w[l] * ek * ek;
				}
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L*N_L) );
}

double dmft_ads_Veff(double f_disp, double D)
{
	double V_eff=0;
	double t[2] = {D, 0};
	
	switch(DOS){
		case 0:
			V_eff = D * D * f_disp * f_disp / 3.0;
			break;
		case 1:
		case 2:
			V_eff = D * D * f_disp * f_disp / 4.0;
			break;
		case 3:
			V_eff = Veff_2d(t);
			break;
		case 4:
			V_eff = Veff_3d(t);
			break;
		case 5:
			V_eff = Veff_4d(t);
			break;
	}
	
	return(V_eff);
}

double dmft_Veff(double D)
{
	return dmft_ads_Veff(1.0, D);
}

static double func_Delta_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = real( ( f_disp*x + V_sqr / (z1 - x - ec) ) / ( z2 - f_disp*x - V_sqr / (z1 - x - ec) ) );
	
	return( f * dos(x, D) );
}
static double func_Delta_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z1 = a[0], z2 = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = imag( ( f_disp*x + V_sqr / (z1 - x - ec) ) / ( z2 - f_disp*x - V_sqr / (z1 - x - ec) ) );
	
	return( f * dos(x, D) );
}

// return : Gc_cav[] ( = Delta[] / V_sqr_imp)
static void Delta_ads(complex<double> *Gc_cav, complex<double> *Gf_perio, complex<double> *self_f, double beta, double D, double V_sqr, double ef, double f_disp, double ec, double chem_pot)
{
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	for(int i=0; i<N_TAU/2; i++){
// 		printf(" %d", i);
		
		complex<double> z1 = IMAG * (double)(2*i+1) * M_PI / beta;
		complex<double> z2 = z1 - ef + chem_pot - self_f[i];
		
		complex<double> a[6] = {z1, z2, V_sqr, f_disp, ec-chem_pot, D};
		
		gsl_function F;
		F.params = a;
		
		double (*func[2])(double, void *) = {func_Delta_re, func_Delta_im};
		double result[2], error[2];
		
// 		for(int j=0; j<201; j++){
// 			double x = -1.0 + (double)j * 0.01;
// 			printf("%d %.5lf %.5lf %.5lf %.5lf\n", j, func[0](x,a), func[1](x,a), func[2](x,a), func[3](x,a));
// 		}
		
		for(int f=0; f<2; f++){
// 			printf("%d", f);
			F.function = func[f];
			
			#if DOS==0 || DOS==1
// 			gsl_integration_qags(&F, -D, D, 0, 1e-7, 1000, w, &result[f], &error[f]);
			gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
			
			#elif DOS==2
			gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
			
			#endif
			
// 			printf("  %.6e %.6e", result[f], error[f]);
		}
// 		printf("\n");
		
		Gc_cav[i] = result[0] + result[1] * IMAG;
	}
	
	gsl_integration_workspace_free(w);
	
	
	double V_sqr_imp = V_sqr + dmft_ads_Veff(f_disp, D);
	
	for(int i=0; i<N_TAU/2; i++){
		Gc_cav[i] /= (Gf_perio[i] * V_sqr_imp);
	}
	
}

//
// input: self_f
//
void dmft_ads1(struct dmft_green_func *G, complex<double> **self_f, int N, double beta, double D, double *V_sqr, double *ef, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c)
{
	
	if( !flag_ensemble ){  // grand canonical ensemble
		
		ave_n = 0;
		for(int s=0; s<N; s++){
			G_site_diag_ads1(G[s].Gf_perio, G[s].Gc_perio, self_f[s], G[s].n_f, G[s].n_c,
			 beta, D, V_sqr[s], ef[s], f_disp, E_c[s], chem_pot);
			ave_n += G[s].n_f + G[s].n_c;
		}
		ave_n /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double chem_pot_old = chem_pot;
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n = 0;
		for(int s=0; s<N; s++){
			G_site_diag_ads1(G[s].Gf_perio, G[s].Gc_perio, self_f[s], G[s].n_f, G[s].n_c,
			 beta, D, V_sqr[s], ef[s], f_disp, E_c[s], chem_pot);
			temp_n += G[s].n_f + G[s].n_c;
		}
		temp_n /= (double)N;
		
		double delta_n = temp_n - ave_n;
		
		int sign_delta = 1;
		if( delta_n < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n) > ACCURACY_FILLING && iter < FILLING_ITER_MAX){
			
			if( (double)sign_delta * delta_n < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n > 0 )  chem_pot -= d_chem_pot;
			else               chem_pot += d_chem_pot;
			
			
			temp_n = 0;
			for(int s=0; s<N; s++){
				G_site_diag_ads1(G[s].Gf_perio, G[s].Gc_perio, self_f[s], G[s].n_f, G[s].n_c,
				 beta, D, V_sqr[s], ef[s], f_disp, E_c[s], chem_pot);
				temp_n += G[s].n_f + G[s].n_c;
			}
			temp_n /= (double)N;
			
			delta_n = temp_n - ave_n;
			
			iter++;
		}
	}
	
	
	for(int s=0; s<N; s++){
		
// 		Delta_ads(G[s].Gc_cav_omega, G[s].Gf_perio, self_f[s],
// 		 beta, D, V_sqr[s], ef[s], f_disp, E_c[s], chem_pot);
		
		double V_sqr_imp = V_sqr[s] + dmft_ads_Veff(f_disp, D);
		
		for(int i=0; i<N_TAU/2; i++){
			complex<double> i_omega = IMAG * (double)(2*i+1) * M_PI / beta;
			
			G[s].Gc_cav_omega[i] = (i_omega - ef[s] + chem_pot - self_f[s][i] - 1.0 / G[s].Gf_perio[i]) / V_sqr_imp;
			
			if(f_disp == 0){
				G[s].self_c[i] = V_sqr[s] / (i_omega - ef[s] + chem_pot - self_f[s][i]);
			}
			else{
				G[s].self_c[i] = 0;
			}
			
			// temp
			G[s].tmat_perio[i] = G[s].Gf_perio[i];
		}
	}
	
}

//============================================================================
//
// for the Anderson lattice model
// dmft_ads2
//   ef -> ef - chem_pot
//   input: Gf
//   f_disp implemented
//

static double func_D_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0], F = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	complex<double> v = V_sqr + (z-x-ec) * x * f_disp;
	double f = real( v / ( (z-x-ec) - v * F ) );
	
	return( f * dos(x, D) );
}
static double func_D_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0], F = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	complex<double> v = V_sqr + (z-x-ec) * x * f_disp;
	double f = imag( v / ( (z-x-ec) - v * F ) );
	
	return( f * dos(x, D) );
}
static double func_Gc2_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0], F = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = real( 1.0 / ( (z-x-ec) - V_sqr * F / (1.0 - f_disp*x*F) ) );
	
	return( f * dos(x, D) );
}
static double func_Gc2_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> z = a[0], F = a[1];
	double V_sqr = real(a[2]), f_disp = real(a[3]), ec = real(a[4]), D = real(a[5]);
	
	double f = imag( 1.0 / ( (z-x-ec) - V_sqr * F / (1.0 - f_disp*x*F) ) );
	
	return( f * dos(x, D) );
}

// input  : G_F[]
// output : G_D[], Gf_perio[], Gc_perio[], n_f, n_c
static void G_site_diag_ads2(complex<double> *Gf_perio, complex<double> *Gc_perio, complex<double> *G_D, complex<double> *G_F, double &n_f, double &n_c, double beta, double D, double V_sqr, double f_disp, double ec, double chem_pot, double jump)
{
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	for(int i=0; i<N_TAU/2; i++){
// 		printf(" %d", i);
		
		complex<double> z = IMAG * (double)(2*i+1) * M_PI / beta;
		
		complex<double> a[6] = {z, G_F[i], V_sqr, f_disp, ec-chem_pot, D};
		
		gsl_function F;
		F.params = a;
		
		double (*func[4])(double, void *) = {func_D_re, func_D_im, func_Gc2_re, func_Gc2_im};
		double result[4], error[4];
		
// 		for(int j=0; j<201; j++){
// 			double x = -1.0 + (double)j * 0.01;
// 			printf("%d %.5lf %.5lf %.5lf %.5lf\n", j, func[0](x,a), func[1](x,a), func[2](x,a), func[3](x,a));
// 		}
		
		for(int f=0; f<4; f++){
// 			printf("%d", f);
			F.function = func[f];
			
			#if DOS==0 || DOS==1
			gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
			
			#elif DOS==2
			gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
			
			#endif
			
// 			printf("  %.6e %.6e", result[f], error[f]);
		}
// 		printf("\n");
		
		
		G_D[i] = result[0] + result[1] * IMAG;
		Gc_perio[i] = result[2] + result[3] * IMAG;
		
		Gf_perio[i] = G_F[i] * (1.0 + G_F[i] * G_D[i]);
	}
	
	gsl_integration_workspace_free(w);
	
	
	// evaluate filling
	n_f = eval_G_tau0m(Gf_perio, beta, jump);
	n_c = eval_G_tau0m(Gc_perio, beta, 1.0);
	
// 	printf(" n_f=%.5lf  n_c=%.5lf  mu=%.5lf\n", n_f, n_c, chem_pot);
	
}

//
// input: Gf
//
/*
void dmft_ads2(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c, double *jump)
{
	
	for(int s=0; s<N; s++){
		double V_sqr_imp = V_sqr[s] + dmft_ads_Veff(f_disp, D);
		
		for(int i=0; i<N_TAU/2; i++){
			G[s].F[i] = Gf[s][i] / (1.0 + Gf[s][i] * G[s].Gc_cav_omega[i] * V_sqr_imp);
		}
	}
	
	complex<double> **G_D;
	array2_alloc(G_D, N, N_TAU/2);
	
	if( !flag_ensemble ){  // grand canonical ensemble
		
		ave_n = 0;
		for(int s=0; s<N; s++){
			G_site_diag_ads2(G[s].Gf_perio, G[s].Gc_perio, G_D[s], G[s].F, G[s].n_f, G[s].n_c,
			 beta, D, V_sqr[s], f_disp, E_c[s], chem_pot, jump[s]);
			ave_n += G[s].n_f + G[s].n_c;
		}
		ave_n /= (double)N;
		
// 		printf("\nFilling : %.5lf\n", n_c);
	}
	else{  // canonical ensemble
		
		double chem_pot_old = chem_pot;
		
		double d_chem_pot = 0.1 * D;
		
		double temp_n = 0;
		for(int s=0; s<N; s++){
			G_site_diag_ads2(G[s].Gf_perio, G[s].Gc_perio, G_D[s], G[s].F, G[s].n_f, G[s].n_c,
			 beta, D, V_sqr[s], f_disp, E_c[s], chem_pot, jump[s]);
			temp_n += G[s].n_f + G[s].n_c;
		}
		temp_n /= (double)N;
		
		double delta_n = temp_n - ave_n;
		
		int sign_delta = 1;
		if( delta_n < 0 )  sign_delta = -1;
		
		int iter=0;
		
		while( fabs(delta_n) > ACCURACY_FILLING && iter < FILLING_ITER_MAX){
			
			if( (double)sign_delta * delta_n < 0 )  sign_delta=0;
			if( sign_delta==0 )  d_chem_pot *= 0.5;
			
			if( delta_n > 0 )  chem_pot -= d_chem_pot;
			else               chem_pot += d_chem_pot;
			
			
			temp_n = 0;
			for(int s=0; s<N; s++){
				G_site_diag_ads2(G[s].Gf_perio, G[s].Gc_perio, G_D[s], G[s].F, G[s].n_f, G[s].n_c,
				 beta, D, V_sqr[s], f_disp, E_c[s], chem_pot, jump[s]);
				temp_n += G[s].n_f + G[s].n_c;
			}
			temp_n /= (double)N;
			
			delta_n = temp_n - ave_n;
			
			iter++;
		}
	}
	
	
	for(int s=0; s<N; s++){
		double V_sqr_imp = V_sqr[s] + dmft_ads_Veff(f_disp, D);
		
		for(int i=0; i<N_TAU/2; i++){
			complex<double> i_omega = IMAG * (double)(2*i+1) * M_PI / beta;
			
			// Delta / V_sqr_eff
			G[s].Gc_cav_omega[i] = G_D[s][i] / (1.0 + G_D[s][i] * G[s].F[i]) / V_sqr_imp;
			
			// temp
// 			G[s].tmat_perio[i] = G[s].Gf_perio[i];
		}
	}
	
	array2_free(G_D, N);
}
*/


struct strct_occup_ads2{
	// input
	complex<double> **G_F;
	int N;
	double beta, D, *V_sqr, f_disp, *ec, chem_pot, *jump, occup_given;
	// output
	complex<double> **Gf_perio, **Gc_perio, **G_D;
	double *n_f, *n_c;
};

static double func_occup_ads2(double d_mu, void *params)
{
	strct_occup_ads2 *p = (strct_occup_ads2 *)params;
	
	for(int s=0; s< p->N; s++){
		G_site_diag_ads2(p->Gf_perio[s], p->Gc_perio[s], p->G_D[s], p->G_F[s], p->n_f[s], p->n_c[s], p->beta, p->D, p->V_sqr[s], p->f_disp, p->ec[s], p->chem_pot + d_mu, p->jump[s]);
	}
	
	double occup_calc = 0;
	for(int s=0; s< p->N; s++)  occup_calc += p->n_f[s] + p->n_c[s];
	occup_calc /= double(p->N);
	
	return occup_calc - p->occup_given;
}

//
// input: Gf
//
void dmft_ads2(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c, double *jump)
{
	for(int s=0; s<N; s++){
		double V_sqr_imp = V_sqr[s] + dmft_ads_Veff(f_disp, D);
		
		for(int i=0; i<N_TAU/2; i++){
			G[s].F[i] = Gf[s][i] / (1.0 + Gf[s][i] * G[s].Gc_cav_omega[i] * V_sqr_imp);
		}
	}
	
	complex<double> **G_D;
	array2_alloc(G_D, N, N_TAU/2);
	
	{
		complex<double> *p_F[N], *p_Gf_perio[N], *p_Gc_perio[N];
		for(int s=0; s<N; s++){
			p_F[s] = G[s].F;
			p_Gf_perio[s] = G[s].Gf_perio;
			p_Gc_perio[s] = G[s].Gc_perio;
		}
		double nf[N], nc[N];
		
		strct_occup_ads2 params = {p_F, N, beta, D, V_sqr, f_disp, E_c, chem_pot, jump, ave_n, p_Gf_perio, p_Gc_perio, G_D, nf, nc};
		
		if( flag_ensemble ){
			chem_pot += tune_chem_pot(func_occup_ads2, (void *)(&params));
		}
		
		func_occup_ads2(0, (void *)(&params));
		for(int s=0; s<N; s++){
			G[s].n_f = nf[s];
			G[s].n_c = nc[s];
		}
	}
	
	for(int s=0; s<N; s++){
		double V_sqr_imp = V_sqr[s] + dmft_ads_Veff(f_disp, D);
		
		for(int i=0; i<N_TAU/2; i++){
			complex<double> i_omega = IMAG * (double)(2*i+1) * M_PI / beta;
			
			// Delta / V_sqr_eff
			G[s].Gc_cav_omega[i] = G_D[s][i] / (1.0 + G_D[s][i] * G[s].F[i]) / V_sqr_imp;
			
			// self_c makes sense only when f_disp=0
			G[s].self_c[i] = G[s].F[i] * V_sqr[s];
		}
	}
	
	array2_free(G_D, N);
}

void dmft_ads2(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c)
{
	double jump[N];
	for(int s=0; s<N; s++)  jump[s] = 1.;
	dmft_ads2(G, Gf, N, beta, D, V_sqr, f_disp, flag_ensemble, ave_n, chem_pot, E_c, jump);
}


//============================================================================
//
// For CPA with hybridization disorder in the periodic Anderson model
//
// dmft_ads3
//   input: tmat, nf, d(nf)/d(mu)
//

struct strct_occup_ads3{
	// input
	complex<double> **self_c;
	int N;
	double beta, D, *ec, chem_pot, occup_given, *nf;
	// output
	complex<double> **Gc_perio;
	double *nc;
};

// n_f fixed
static double func_occup_ads3(double d_mu, void *params)
{
	strct_occup_ads3 *p = (strct_occup_ads3 *)params;
	
	for(int s=0; s< p->N; s++){
// 		G[a].n_c = Gc_site_diag(G[a].Gc_perio, G[a].self_c, beta, chem_pot-E_c[a], D);
		p->nc[s] = Gc_site_diag(p->Gc_perio[s], p->self_c[s], p->beta, p->chem_pot + d_mu - p->ec[s], p->D);
	}
	
	double occup_calc = 0;
	for(int s=0; s< p->N; s++)  occup_calc += p->nf[s] + p->nc[s];
	occup_calc /= double(p->N);
	
	return occup_calc - p->occup_given;
}

// for CPA with hybridization disorder
//
// input: tmat, nf
//
void dmft_ads3(struct dmft_green_func *G, complex<double> **tmat, double *nf_in, int N, double beta, double D, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c)
{
	for(int s=0; s<N; s++){
		for(int i=0; i<N_TAU/2; i++){
			G[s].self_c[i] = tmat[s][i] / (1.0 + tmat[s][i] * G[s].Gc_cav_omega[i]);
		}
	}
	
	{
		complex<double> *p_self_c[N], *p_Gc_perio[N];
		double nc[N], nf[N];
		for(int s=0; s<N; s++){
			p_self_c[s] = G[s].self_c;
			p_Gc_perio[s] = G[s].Gc_perio;
			nf[s] = nf_in[s];
		}
		
		strct_occup_ads3 params = {p_self_c, N, beta, D, E_c, chem_pot, ave_n, nf, p_Gc_perio, nc};
		
		if( flag_ensemble ){
			chem_pot += tune_chem_pot(func_occup_ads3, (void *)(&params));
		}
		func_occup_ads3(0, (void *)(&params));
		
		for(int s=0; s<N; s++){
			G[s].n_f = nf[s];
			G[s].n_c = nc[s];
		}
		if( !flag_ensemble ){
			ave_n = 0;
			for(int s=0; s<N; s++)  ave_n += (nf[s] + nc[s]) / double(N);
		}
	}
	
	for(int s=0; s<N; s++){
		for(int i=0; i<N_TAU/2; i++){
			G[s].Gc_cav_omega[i] = G[s].Gc_perio[i] / ( 1. + G[s].Gc_perio[i] * G[s].self_c[i] );
			G[s].tmat_perio[i] = G[s].self_c[i] * ( 1. + G[s].self_c[i] * G[s].Gc_perio[i] );
		}
	}
}


//============================================================================
//
// Hubbard model
//

/*
static double func_hubbard_re(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> F_inv = a[0];
	double delta_mu = real(a[1]), D = real(a[2]);
	
	double f = real( 1./ ( F_inv + delta_mu - x ) );
	
	return( f * dos(x, D) );
}
static double func_hubbard_im(double x, void *params)
{
	complex<double> *a = (complex<double> *)params;
	complex<double> F_inv = a[0];
	double delta_mu = real(a[1]), D = real(a[2]);
	
	double f = imag( 1./ ( F_inv + delta_mu - x ) );
	
	return( f * dos(x, D) );
}

// input  : G_F[]
// output : Gf_perio[], n_f
static void G_site_diag_hubbard(complex<double> *Gf_perio, complex<double> *G_F, double &n_f, double beta, double D, double delta_mu, double jump)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	for(int i=0; i<N_TAU/2; i++){
// 		printf(" %d", i);
		
// 		complex<double> z = IMAG * (double)(2*i+1) * M_PI / beta;
		
		complex<double> a[3] = {1./ G_F[i], delta_mu, D};
		
		gsl_function F;
		F.params = a;
		
		double (*func[2])(double, void *) = {func_hubbard_re, func_hubbard_im};
		double result[2], error[2];
		
// 		for(int j=0; j<201; j++){
// 			double x = -1.0 + (double)j * 0.01;
// 			printf("%d %.5lf %.5lf %.5lf %.5lf\n", j, func[0](x,a), func[1](x,a), func[2](x,a), func[3](x,a));
// 		}
		
		for(int f=0; f<2; f++){
// 			printf("%d", f);
			F.function = func[f];
			
			#if DOS==0 || DOS==1
			gsl_integration_qag(&F, -D, D, GSL_EPSABS, GSL_EPSREL, 1000, GSL_KEY, w, &result[f], &error[f]);
			
			#elif DOS==2
			gsl_integration_qagi(&F, GSL_EPSABS, GSL_EPSREL, 1000, w, &result[f], &error[f]);
			
			#endif
			
// 			printf("  %.6e %.6e", result[f], error[f]);
		}
// 		printf("\n");
		
		
		Gf_perio[i] = result[0] + result[1] * IMAG;
	}
	
	gsl_integration_workspace_free(w);
	
	
	// evaluate filling
	n_f = eval_G_tau0m(Gf_perio, beta, jump);
// 	double Gf_tau_temp[N_TAU];
// 	fft_fermion_radix2_omega2tau(Gf_tau_temp, Gf_perio, beta, N_TAU);
// 	n_f = 1.0 + Gf_tau_temp[0];
}


static void G_site_diag_hubbard(complex<double> *Gf_perio, complex<double> *G_F, double &n_f, double beta, double D, double delta_mu, double jump)
{
	for(int i=0; i<N_TAU/2; i++){
		complex<double> z = 1./ G_F[i] + delta_mu;
		
		Gf_perio[i] = hilbert_transform(z, D);
	}
	
	// evaluate filling
	n_f = eval_G_tau0m(Gf_perio, beta, jump);
// 	double Gf_tau_temp[N_TAU];
// 	fft_fermion_radix2_omega2tau(Gf_tau_temp, Gf_perio, beta, N_TAU);
// 	n_f = 1.0 + Gf_tau_temp[0];
}
*/

struct params_dmft_hubbard{
	struct dmft_green_func *G;
	double beta;
	double D;
	int N;
	double occup;
// 	double chem_pot;
	double *jump;
	complex<double> **Gf;  // used only when DOS==1 (Bethe lattice)
};

// input: delta_mu
// return: deviation of occupation
static double func_dmft_hubbard(double delta_mu, void *params)
{
	struct params_dmft_hubbard *p = (params_dmft_hubbard *)params;
	
// 	double r = - p->occup;
// 	for(int s=0; s<p->N; s++){
// 		G_site_diag_hubbard(p->G[s].Gf_perio, p->G[s].F, p->G[s].n_f, p->beta, p->D, delta_mu);
// 		r += p->G[s].n_f;
// 	}
// 	
// 	return r;
	
// 	double occup = 0;
// 	for(int s=0; s<p->N; s++){
// 		G_site_diag_hubbard(p->G[s].Gf_perio, p->G[s].F, p->G[s].n_f, p->beta, p->D, delta_mu, p->jump[s]);
// 		occup += p->G[s].n_f;
// 	}
// 	if( p->N == 1 )  occup *= 2.;
// 	
// 	return( occup - p->occup );
	
	double occup = 0;
	for(int s=0; s<p->N; s++){
		
		// Hilbert transform
		for(int i=0; i<N_TAU/2; i++){
			complex<double> z = 1./ p->G[s].F[i] + delta_mu;
			p->G[s].Gf_perio[i] = hilbert_transform(z, p->D);
		}
		
		// evaluate filling
		p->G[s].n_f = eval_G_tau0m(p->G[s].Gf_perio, p->beta, p->jump[s]);
		occup += p->G[s].n_f;
	}
	
	if( DOS == 1 ){  // Bethe lattice
		occup = 0;
		for(int s=0; s<p->N; s++){
			complex<double> Gf_new[N_TAU/2];
			for(int i=0; i<N_TAU/2; i++){
				Gf_new[i] = p->Gf[s][i] / (1. + p->Gf[s][i] * delta_mu);
			}
			occup += eval_G_tau0m(Gf_new, p->beta, p->jump[s]);
		}
	}
	
	if( p->N == 1 )  occup *= 2.;
	
	return( occup - p->occup );
}
//
// input: Gf
//
void dmft_hubbard(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, int flag_ensemble, double &occup, double &chem_pot, double *jump)
{
	double V_sqr = dmft_Veff(D);
	
	for(int s=0; s<N; s++){
		for(int i=0; i<N_TAU/2; i++){
			G[s].F[i] = Gf[s][i] / (1.0 + Gf[s][i] * G[s].Gc_cav_omega[i] * V_sqr);
		}
	}
	
	double delta_mu = 0;
	
	if( !flag_ensemble ){  // chemical potential given
		struct params_dmft_hubbard params = {G, beta, D, N, 0, jump, Gf};
		occup = func_dmft_hubbard(delta_mu, &params);
	}
	else{  // occupation given
		struct params_dmft_hubbard params = {G, beta, D, N, occup, jump, Gf};
		
		gsl_function F;
		F.function = &func_dmft_hubbard;
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
			delta_mu = gsl_root_fsolver_root (s);
			x_lo = gsl_root_fsolver_x_lower (s);
			x_hi = gsl_root_fsolver_x_upper (s);
			status = gsl_root_test_interval (x_lo, x_hi, 0, ACCURACY_FILLING);
			
// 			if (status == GSL_SUCCESS)  printf ("Converged:\n");
// 			printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, delta_mu, x_hi - x_lo);
		}
		while (status == GSL_CONTINUE && iter < FILLING_ITER_MAX);
		
// 		printf("%s", status);
		gsl_root_fsolver_free (s);
		
		chem_pot += delta_mu;
	}
	
	if( DOS == 1 ){  // Bethe lattice
		for(int s=0; s<N; s++){
			for(int i=0; i<N_TAU/2; i++){
				// Delta / V_sqr
// 				G[s].Gc_cav_omega[i] = 0.25 * D * D * G[s].Gf_perio[i] / V_sqr;
				
				G[s].Gc_cav_omega[i] = 0.25 * D * D * Gf[s][i] / V_sqr;
				
// 				complex<double> Gf_new = Gf[s][i] / (1. + Gf[s][i] * delta_mu);
// 				G[s].Gc_cav_omega[i] = 0.25 * D * D * Gf_new / V_sqr;
			}
		}
	}
	else{
		for(int s=0; s<N; s++){
			for(int i=0; i<N_TAU/2; i++){
				// Delta / V_sqr
				G[s].Gc_cav_omega[i] = ( 1./ G[s].F[i] - 1./ G[s].Gf_perio[i] + delta_mu) / V_sqr;
			}
		}
	}
}
void dmft_hubbard(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, int flag_ensemble, double &occup, double &chem_pot)
{
	double jump[N];
	for(int s=0; s<N; s++)  jump[s] = 1;
	dmft_hubbard(G, Gf, N, beta, D, flag_ensemble, occup, chem_pot, jump);
}

// layer Hubbard model
void dmft_layer_hb(struct dmft_green_func **G, complex<double> **Gf, int N, double beta, double D, const CPPL::dsymatrix &hop)
// void dmft_layer_hb(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double **hop)
{
	double V_sqr = dmft_Veff(D);
	
	for(int l=0; l<N; l++){
		for(int i=0; i<N_TAU/2; i++){
			G[l]->F[i] = Gf[l][i] / (1.0 + Gf[l][i] * G[l]->Gc_cav_omega[i] * V_sqr);
		}
	}
	
	CPPL::zgematrix mat(N,N), ur, ur_inv, Glat;
	vector< complex<double> > xi, g_bar(N);
	for(int i=0; i<N_TAU/2; i++){
		mat.zero();
		mat += hop.to_zhematrix();
		for(int l=0; l<N; l++)  mat(l,l) += 1./ G[l]->F[i];
		
		CPPL_zgeev(mat, xi, ur, ur_inv);
		
		// Hilbert transform
		for(int l=0; l<N; l++)
			g_bar[l] = hilbert_transform(xi[l], D);
		
		Glat = ur * g_bar * ur_inv;
		
		for(int l=0; l<N; l++)  G[l]->Gf_perio[i] = Glat(l,l);
	}
	
	// evaluate filling
	for(int l=0; l<N; l++){
		G[l]->n_f = eval_G_tau0m(G[l]->Gf_perio, beta);
	}
	
	
// 	if( DOS == 1 ){  // Bethe lattice
// 		for(int s=0; s<N; s++){
// 			for(int i=0; i<N_TAU/2; i++){
// 				G[s].Gc_cav_omega[i] = 0.25 * D * D * Gf[s][i] / V_sqr;
// 			}
// 		}
// 	}
// 	else{
		for(int l=0; l<N; l++){
			for(int i=0; i<N_TAU/2; i++){
				// Delta / V_sqr
				G[l]->Gc_cav_omega[i] = ( 1./ G[l]->F[i] - 1./ G[l]->Gf_perio[i]) / V_sqr;
// 				G[l]->Gc_cav_omega[i] = ( 1./ G[l]->F[i] - 1./ G[l]->Gf_perio[i] + delta_mu) / V_sqr;
			}
		}
// 	}
}

//============================================================================
//
// Two-particle quantities
//

// Evaluate P_{12}(w,nu) = - (1/N) \sum_k g_1(k,w) g_2(k,w+nu)
// input:
//   xi1, xi2,
//   g1 = hilbert_transform(xi1, D);
//   g2 = hilbert_transform(xi2, D);
static complex<double> two_particle_unif(complex<double> xi1, complex<double> xi2, complex<double> g1, complex<double> g2, double D)
{
	if( abs(xi1 - xi2) > 1e-8 ){  // xi1 != xi2
		return (g1 - g2) / (xi1 - xi2);
	}
	else{
		return hilbert_sqr_transform(xi1, D);
	}
}

// Evaluate the corrent-corrent correlation function K[nu]
// input: G[l]->F[i], which has been computed in dmft_layer_hb()
void dmft_layer_current(struct dmft_green_func **G, int N, double beta, double D, const CPPL::dsymatrix &hop, double *Kz, double **Kxy, int N_nu)
{
	int nw = N_TAU/2;
	
	CPPL::zgematrix mat(N,N);
	CPPL::zgematrix ur[nw], ur_inv[nw];
	vector< complex<double> > xi[nw], g_bar[nw];
	for(int i=0; i<nw; i++){  // MPI
		mat.zero();
		mat += hop.to_zhematrix();
		for(int l=0; l<N; l++)  mat(l,l) += 1./ G[l]->F[i];
		
		CPPL_zgeev(mat, xi[i], ur[i], ur_inv[i]);
		
		// Hilbert transform
		g_bar[i].resize(N);
		for(int l=0; l<N; l++)
			g_bar[i][l] = hilbert_transform(xi[i][l], D);
	}
	
	// P_{c1,c2}(w,nu) = - (1/N) \sum_{k} g_{c1}(k,w) g_{c2}(k,w+nu)
	// P[nu][w](c1,c2)
	CPPL::zgematrix P[N_nu][nw];
	
	for(int nu=0; nu<N_nu; nu++){
		for(int i=0; i<nw-nu; i++){  // MPI
			P[nu][i].resize(N,N);
			for(int c1=0; c1<N; c1++){
				for(int c2=0; c2<N; c2++){
					P[nu][i](c1,c2) = 
					  two_particle_unif(xi[i][c1], xi[i+nu][c2], g_bar[i][c1], g_bar[i+nu][c2], D);
				}
			}
		}
		
		// current in the z direction
		Kz[nu]=0;
		for(int a1=0; a1<N; a1++){
			for(int b1=a1+1; b1<N; b1++){  // b1 > a1
				if( hop(b1,a1) != 0 ){
					for(int a2=0; a2<N; a2++){
						for(int b2=a2+1; b2<N; b2++){  // b2 > a2
							if( hop(b2,a2) != 0 ){  // hop(b2,a2) == hop(a2,b2) is assumed
								complex<double> temp = 0;
								for(int i=0; i<nw-nu; i++){  // MPI
									for(int c1=0; c1<N; c1++){
										for(int c2=0; c2<N; c2++){
											// Pi_{a1,b1;b2,a2} - Pi_{a1,b1;a2,b2}
											// Pi = ur[i]() * ur_inv[i]() * p() * ur[i+nu]() * ur_inv[i+nu]()
											temp += ( ur[i](a2,c1) * ur_inv[i+nu](c2,b2) - ur[i](b2,c1) * ur_inv[i+nu](c2,a2) )
											       * ur_inv[i](c1,a1) * ur[i+nu](b1,c2) * P[nu][i](c1,c2);
										}
									}
								}
								Kz[nu] += hop(b1,a1) * hop(b2,a2) * real(temp);
							}
						}
					}
				}
			}
		}
		Kz[nu] *= 2. / beta / double(N);
		
		// current in the layer
		for(int a=0; a<N; a++){
			complex<double> temp = 0;
			for(int i=0; i<nw-nu; i++){  // MPI
				for(int c1=0; c1<N; c1++){
					for(int c2=0; c2<N; c2++){
						// Pi_{a,a;a,a}
						temp += ur[i](a,c1) * ur_inv[i+nu](c2,a)
						       * ur_inv[i](c1,a) * ur[i+nu](a,c2) * P[nu][i](c1,c2);
					}
				}
			}
			Kxy[a][nu] = (2. / beta) * real(temp);
			Kxy[a][nu] *= D * D;  // adjust dimension (K: energy)
		}
	}
}

void dmft_two_particle(struct dmft_green_func &G, struct dmft_tp &T0, double beta, double D, int n_f, double chem_pot)
{
	
	
	#if DOS==2 // hyper-cubic
	
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		complex<double> z = IMAG * omega_n - G.self_c[i] + chem_pot;
		
		T0.Pi_loc[i] = -G.Gc_perio[i] * G.Gc_perio[i];
		T0.Pi_unif[i] = -4.0 / D / D * (z * G.Gc_perio[i] - 1.0);
		T0.Pi_stag[i] = -G.Gc_perio[i] / z;
	
	}
	
	// static susceptibility
	T0.stat_Pi_loc = 0;
	T0.stat_Pi_unif = 0;
	T0.stat_Pi_stag = 0;
	
	for(int i=0; i<N_TAU/2; i++){
		double suscep0 = 1.0 / pow((double)(2*i+1) * M_PI / beta, 2);
		
		T0.stat_Pi_loc += real(T0.Pi_loc[i]) - suscep0;
		T0.stat_Pi_unif += real(T0.Pi_unif[i]) - suscep0;
		T0.stat_Pi_stag += real(T0.Pi_stag[i]) - suscep0;
	}
	
	T0.stat_Pi_loc *= 2.0 / beta;
	T0.stat_Pi_unif *= 2.0 / beta;
	T0.stat_Pi_stag *= 2.0 / beta;
	
	T0.stat_Pi_loc += 0.25 * beta;
	T0.stat_Pi_unif += 0.25 * beta;
	T0.stat_Pi_stag += 0.25 * beta;
	
	
	#endif
	
	for(int i=0; i<N_TAU/2; i++){
		complex<double> self_sqr = G.self_c[i] * G.self_c[i];
		complex<double> self_4power = self_sqr * self_sqr;
		
		complex<double> temp = -self_sqr * ( 1.0 + 2.0 * G.self_c[i] * G.Gc_perio[i] );
		
		T0.T0_loc[i] = (double)n_f * (temp + self_4power * T0.Pi_loc[i]);
		T0.T0_unif[i] = (double)n_f * (temp + self_4power * T0.Pi_unif[i]);
		T0.T0_stag[i] = (double)n_f * (temp + self_4power * T0.Pi_stag[i]);
	}
	
}

void dmft_two_particle_dynamic(struct dmft_green_func &G, complex<double> *Pi_unif_dynamic, double beta, double D, double chem_pot)
{
	
	complex<double> Gc[N_TAU], self[N_TAU];
	
	for(int i=0; i<N_TAU/2; i++){
		Gc[i] = G.Gc_perio[i];
		self[i] = G.self_c[i];
	}
	for(int i=N_TAU/2; i<N_TAU; i++){
		Gc[i] = 1.0 / (IMAG * (double)(2*i+1) * M_PI / beta);
		self[i] = 1.0 / (IMAG * (double)(2*i+1) * M_PI / beta);
	}
	
// 	double Pi_unif_dynamic[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Pi_unif_dynamic[i] = 0;
	
	
	// i == 0
	
	for(int j=0; j<N_TAU/2; j++){
		double omega_n = (double)(2*j+1) * M_PI / beta;
		complex<double> z = IMAG * omega_n - G.self_c[j] + chem_pot;
		
		Pi_unif_dynamic[0] += real( -4.0 / D / D * (z * G.Gc_perio[j] - 1.0) ) - 1.0 / pow(omega_n, 2);
	}
	
	
	// i != 0
	
	for(int i=1; i<N_TAU/2; i++){
		
		complex<double> i_nu = IMAG * (double)(2*i) * M_PI / beta;
		
		int i2 = i/2;
		
		if(i%2==1){
			Pi_unif_dynamic[i] += 0.5 * real((conj(Gc[i2]) - Gc[i2]) / (-i_nu - conj(self[i2]) + self[i2]))
			                    - 0.5 * pow(beta/M_PI, 2) / (double)(-(2*i2+1) * (2*i2+1));
			
// 			printf("%d %.3e\n", i, 0.5 * real((conj(Gc[i2]) - Gc[i2]) / (-i_nu - conj(self[i2]) + self[i2])) );
		}
		
		for(int j=0; j<i2; j++){
			int j2 = i - 1 - j;
			
			Pi_unif_dynamic[i] += real( (conj(Gc[j]) - Gc[j2]) / (-i_nu - conj(self[j]) + self[j2]) )
			              - pow(beta/M_PI, 2) / (double)(-(2*j+1) * (2*j2+1));
		}
		for(int j=0; j<N_TAU/2; j++){
			int j2 = j + i;
			
			Pi_unif_dynamic[i] += real( (Gc[j] - Gc[j2]) / (-i_nu - self[j] + self[j2]) )
			              - pow(beta/M_PI, 2) / (double)((2*j+1) * (2*j2+1));
		}
	}
	
	for(int i=0; i<N_TAU/2; i++){
		Pi_unif_dynamic[i] *= 2.0 / beta;
	}
	
	Pi_unif_dynamic[0] += 0.25 * beta;
	
}


//============================================================================

static void print_dmft_sub1(struct dmft_green_func *G, int N, double pot_scat, int num, int subl, double mu)
{
	char filename[256];
	if(subl<0)  sprintf(filename, DATA_DIR "%02d-dmft.dat", num);
	else        sprintf(filename, DATA_DIR "%02d-%c-dmft.dat", num, subl+'A');
	
	FILE *fp=fopen(filename, "w");
	if( ! isnan(mu) )  fprintf(fp, "# mu=%.6e\n", mu);
	
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		
		for(int a=0; a<N; a++){
			complex<double> Gc_cav = G[a].Gc_cav_omega[i] / ( 1.0 + pot_scat * G[a].Gc_cav_omega[i] );
			
			complex<double> tmat = pot_scat * (1.0 + pot_scat * G[a].Gc_cav_omega[i])
			 + G[a].tmat_perio[i] * pow(1.0 + pot_scat * G[a].Gc_cav_omega[i], 2);
			
			
			fprintf(fp, "  %.8e %.8e %.5e %.5e %.5e %.5e %.5e %.5e",
			 real(Gc_cav), imag(Gc_cav),
			 real(G[a].Gc_perio[i]), imag(G[a].Gc_perio[i]),
			 real(G[a].self_c[i]+pot_scat), imag(G[a].self_c[i]),
			 real(tmat), imag(tmat));
		}
		
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	printf("\n '%s'\n", filename);
}
static void print_dmft_sub2(struct dmft_green_func *G, int N, double pot_scat, char *filename)
{
	FILE *fp=fopen(filename, "w");
	
// 	for(int i=0; i<N_TAU/2; i++){
	for(int i=0; i<30; i++){
		
		fprintf(fp, "%d", i);
		
		for(int a=0; a<N; a++){
			complex<double> Gc_cav = G[a].Gc_cav_omega[i] / ( 1.0 + pot_scat * G[a].Gc_cav_omega[i] );
			
// 			fprintf(fp, " %.3e %.3e %.3e %.3e",
// 			 real(Gc_cav), imag(Gc_cav), real(G[a].self_c[i]+pot_scat), imag(G[a].self_c[i]));
			
			fprintf(fp, " %.3e %.3e", real(Gc_cav), imag(Gc_cav));
		}
		
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	printf("\n '%s'\n", filename);
}


void print_dmft_num(struct dmft_green_func *G, int N, int num1, int num2)
{
	print_dmft_num(G, N, 0.0, num1, num2);
}
void print_dmft_num(struct dmft_green_func *G, int N, double pot_scat, int num1, int num2)
{
	char filename[256];
	sprintf(filename, DATA_DIR "%02d-dmft-%03d.dat", num1, num2);
	
	print_dmft_sub2(G, N, pot_scat, filename);
}


void print_dmft(struct dmft_green_func *G, int N, double pot_scat, int num, int subl)
{
	print_dmft_sub1(G, N, pot_scat, num, subl, NAN);
}
void print_dmft(struct dmft_green_func *G, int N, double pot_scat, int num)
{
	print_dmft_sub1(G, N, pot_scat, num, -1, NAN);
}
void print_dmft(struct dmft_green_func *G, int N, int num, int subl)
{
	print_dmft_sub1(G, N, 0.0, num, subl, NAN);
}
void print_dmft(struct dmft_green_func *G, int N, int num)
{
	print_dmft_sub1(G, N, 0.0, num, -1, NAN);
}

void print_dmft_mu(struct dmft_green_func *G, int N, int num, int subl, double mu)
{
	print_dmft_sub1(G, N, 0.0, num, subl, mu);
}
void print_dmft_mu(struct dmft_green_func *G, int N, int num, double mu)
{
	print_dmft_sub1(G, N, 0.0, num, -1, mu);
}


// 1: G is set
// 0: G is not set
static int read_dmft_sub(struct dmft_green_func *G, int N, int num, int subl, double &mu)
{
	char filename[256];
	if(subl<0)  sprintf(filename, "%02d-dmft.dat", num);
	else        sprintf(filename, "%02d-%c-dmft.dat", num, subl+'A');
	
	FILE *fp=fopen(filename, "r");
	if( fp == NULL )  return 0;
	
	{  // chemical potential
		char str[1024];
		if( fgets(str, 1024, fp) == NULL)  return 0;
		if( strncmp(str, "# mu=", 5) == 0 )  sscanf(str, "# mu=%lf", &mu);
// 		printf(" mu = %.6e\n", mu);
		rewind(fp);  // move to the head of the file
	}
	
	{
// 		if(my_rank==0){
// 			printf("\n---\n");
// 			printf("\n '%s' opened\n", filename);
// 		}
		
		int n=0;  // data number
		double data[2048];
		
		while( read_data(fp, data) >= 1+8*N ){
			for(int s=0; s<N; s++){
				G[s].Gc_cav_omega[n] = data[1+8*s] + IMAG * data[2+8*s];
				G[s].Gc_perio[n]     = data[3+8*s] + IMAG * data[4+8*s];
				G[s].self_c[n]       = data[5+8*s] + IMAG * data[6+8*s];
				G[s].tmat_perio[n]   = data[7+8*s] + IMAG * data[8+8*s];
			}
			n++;
		}
		fclose(fp);
		
		// test
// 		for(int i=0; i<n; i++){
// 			printf("%d", i);
// 			for(int a=0; a<N_F; a++)  printf(" %.5e %.5e", real(G0_omega[a][i]), imag(G0_omega[a][i]));
// 			printf("\n");
// 		}
		
		
		if(n != N_TAU/2){
// 			if(my_rank==0){
				printf("\n*** error : inconsistency in data number\n");
				printf("  n = %d,  N_TAU/2 = %d\n", n, N_TAU/2);
// 			}
			exit(0);
		}
// 		if(my_rank==0)  printf(" G0_omega[] has been set\n");
	}
	return 1;
}

int read_dmft(struct dmft_green_func *G, int N, int num, int subl)
{
	double dummy;
	return read_dmft_sub(G, N, num, subl, dummy);
}
int read_dmft(struct dmft_green_func *G, int N, int num)
{
	double dummy;
	return read_dmft_sub(G, N, num, -1, dummy);
}
int read_dmft_mu(struct dmft_green_func *G, int N, int num, int subl, double &mu)
{
	return read_dmft_sub(G, N, num, subl, mu);
}
int read_dmft_mu(struct dmft_green_func *G, int N, int num, double &mu)
{
	return read_dmft_sub(G, N, num, -1, mu);
}



void print_dmft_tp(struct dmft_tp &T0, int num)
{
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_tp.dat", num);
	
	FILE *fp=fopen(filename, "w");
	
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		
		fprintf(fp, " %.6e %.6e", real(T0.Pi_loc[i]), imag(T0.Pi_loc[i]));
		fprintf(fp, " %.6e %.6e", real(T0.Pi_unif[i]), imag(T0.Pi_unif[i]));
		fprintf(fp, " %.6e %.6e", real(T0.Pi_stag[i]), imag(T0.Pi_stag[i]));
		
		fprintf(fp, " %.6e %.6e", real(T0.T0_loc[i]), imag(T0.T0_loc[i]));
		fprintf(fp, " %.6e %.6e", real(T0.T0_unif[i]), imag(T0.T0_unif[i]));
		fprintf(fp, " %.6e %.6e", real(T0.T0_stag[i]), imag(T0.T0_stag[i]));
		fprintf(fp, " %.6e %.6e", real(T0.T0_loc2[i]), imag(T0.T0_loc2[i]));
		
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}

void print_dmft_tp_dynamic(complex<double> *Pi_unif_dynamic, double beta, double D, int num)
{
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_tp_dynamic.dat", num);
	
	FILE *fp=fopen(filename, "w");
	
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		
		fprintf(fp, " %.5e", real(Pi_unif_dynamic[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	complex<double> Pi_pade[N_TAU/2], i_omega_b[N_TAU/2];
	
	for(int i=0; i<N_TAU/2; i++)  i_omega_b[i] = IMAG * (double)(2*i) * M_PI / beta;
	pade_init(i_omega_b, Pi_unif_dynamic, Pi_pade, N_TAU/2);
	
	sprintf(filename, DATA_DIR "%02d-dmft_tp_dynamic_pade.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<1000; i++){
		complex<double> w = 1.0 * D / 999. * (double)i;
		complex<double> temp = pade_calc(i_omega_b, Pi_pade, w, N_TAU/2);
		
		fprintf(fp, "%.5e %.5e %.5e", real(w), real(temp), imag(temp));
		
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}


static double dmft_nk_cs_sub1(complex<double> *z, double ek, double beta)
{
	// occupation number n(e)
	
	complex<double> Gk[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Gk[i] = 1.0 / (z[i] - ek);
	
	double sum_Gk=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk += real(Gk[i]);
	}
	sum_Gk *= 2.0 / beta;
	sum_Gk += 0.5;
	
	return( sum_Gk );
}
static double dmft_nk_cs_sub2(complex<double> *z, double ek, double beta)
{
	// -dn(e)/de
	
	complex<double> Gk[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Gk[i] = 1.0 / (z[i] - ek);
	
	double sum_Gk_sqr=0;
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk_sqr += real(Gk[i] * Gk[i]) + 1.0 / (omega_n*omega_n);
	}
	sum_Gk_sqr *= 2.0 / beta;
	sum_Gk_sqr -= 0.25 * beta;
	
	return( -sum_Gk_sqr );
}
void dmft_nk_cs(struct dmft_green_func &G, double beta, double D, double chem_pot, double nc, int num)
{
	complex<double> z[N_TAU/2];
	
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		z[i] = IMAG * omega_n - G.self_c[i] + chem_pot;
	}
	
	
	int N_p=1001;
	double ek[N_p], ne[N_p], ne_de[N_p];
	
	double ek_max=1.5*D;
	for(int ik=0; ik<N_p; ik++){
		ek[ik] = -ek_max + 2.0 * ek_max * (double)ik / (double)(N_p-1);
		
		ne[ik] = dmft_nk_cs_sub1(z, ek[ik], beta);
		ne_de[ik] = dmft_nk_cs_sub2(z, ek[ik], beta);
	}
	
	
	double mu_small = eval_chem_pot(nc, D);
	double mu_large = eval_chem_pot(nc + 0.5, D);
	
	double ne_small = dmft_nk_cs_sub1(z, mu_small, beta);
	double ne_large = dmft_nk_cs_sub1(z, mu_large, beta);
	double ne_de_small = dmft_nk_cs_sub2(z, mu_small, beta);
	double ne_de_large = dmft_nk_cs_sub2(z, mu_large, beta);
	
	//
	// find local maxima of ne_de
	//  using the golden section search
	//
	double ek_peak[3], ne_peak[3], ne_de_peak[3];
	int n_peak=0;
	
	int flag_pm=+1;
	
	for(int ik=1; ik<N_p; ik++){
		if( ne_de[ik+1] < ne_de[ik] ){
			
			if(flag_pm){
// 				printf(" %.5e %.5e %.5e\n", ek[ik], ne[ik], ne_de[ik]);
				
				double x1=ek[ik-1], x2=ek[ik], x3=ek[ik+1];
				double y1=ne_de[ik-1], y2=ne_de[ik], y3=ne_de[ik];
				
				double golden = 0.5*(3.0-sqrt(5.0));
				
				do{
					if( x2-x1 > x3-x2 ){
						double x4 = x2 - golden * (x2-x1);
						double y4 = dmft_nk_cs_sub2(z, x4, beta);
						
						if(y4<y2){
							x1 = x4;
							y1 = y4;
						}
						else{
							x3 = x2;
							y3 = y2;
							x2 = x4;
							y2 = y4;
						}
						
					}
					else{
						double x4 = x2 + golden * (x3-x2);
						double y4 = dmft_nk_cs_sub2(z, x4, beta);
						
						if(y4<y2){
							x3 = x4;
							y3 = y4;
						}
						else{
							x1 = x2;
							y1 = y2;
							x2 = x4;
							y2 = y4;
						}
					}
					
				}while( x3-x1 > 1.0e-5 );
				
				ek_peak[n_peak] = x2;
				ne_peak[n_peak] = dmft_nk_cs_sub1(z, x2, beta);
				ne_de_peak[n_peak] = y2;
				
				n_peak++;
				if(n_peak==3) break;
			}
			
			flag_pm = 0;
		}
		else{
			flag_pm = 1;
		}
	}
	
	
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_ne.dat", num);
	
	FILE *fp=fopen(filename, "w");
	fprintf(fp, "# small  %.6e %.6e %.6e\n", mu_small, ne_small, ne_de_small);
	fprintf(fp, "# large  %.6e %.6e %.6e\n", mu_large, ne_large, ne_de_large);
	
	for(int ik=0; ik<n_peak; ik++){
		fprintf(fp, "# peak%d  %.6e %.6e %.6e\n", ik, ek_peak[ik], ne_peak[ik], ne_de_peak[ik]);
	}
	for(int ik=0; ik<N_p; ik++){
		fprintf(fp, "%.6e %.6e %.6e\n", ek[ik], ne[ik], ne_de[ik]);
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	
// 	#if DOS==2
// 	
// 	int N_nu=5;
// 	double p[N_p], nk[N_nu][N_p], nk_dk[N_nu][N_p];
// 	
// 	for(int nu=0; nu<N_nu; nu++){
// 		
// 		double ep[N_p];
// 		for(int ik=0; ik<N_p; ik++){
// 			p[ik] = (double)ik / (double)(N_p-1);
// 			ep[ik] = -D / sqrt(2.0) * cos(p[ik] * M_PI) * (double)(1+nu);
// 			
// 			nk[nu][ik] = dmft_nk_cs_sub1(z, ek[ik], beta);
// 			nk_dk[nu][ik] = dmft_nk_cs_sub2(z, ek[ik], beta);
// 			
// 			nk_dk[nu][ik] *= D / sqrt(2.0) * sin(p[ik] * M_PI) * (double)(1+nu);
// 		}
// 	}
// 	
// 	
// 	sprintf(filename, DATA_DIR "%02d-dmft_nk.dat", num);
// 	
// 	fp=fopen(filename, "w");
// 	for(int ik=0; ik<N_p; ik++){
// 		fprintf(fp, "%.6e", p[ik]);
// 		for(int nu=0; nu<N_nu; nu++)  fprintf(fp, " %.6e %.6e", nk[nu][ik], nk_dk[nu][ik]);
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
// 	
// 	#endif
	
}

void dmft_nk_cs(struct dmft_green_func *G, int N, double beta, double D, double chem_pot, double nc, int num)
{
	double ec[N];
	for(int a=0; a<N; a++)  ec[a] = 0;
	
	dmft_nk_cs(G, N, beta, D, chem_pot, nc, ec, num);
}

void dmft_nk_cs(struct dmft_green_func *G, int N, double beta, double D, double chem_pot, double nc, double *ec, int num)
{
	int N_p=1001;
	double ek[N_p], ne[N][N_p], ne_de[N][N_p];
	
	double mu_small[N], mu_large[N];
	double ne_small[N], ne_large[N], ne_de_small[N], ne_de_large[N];
	
	double ek_peak[N][3], ne_peak[N][3], ne_de_peak[N][3];
	int n_peak[N];
	
	for(int a=0; a<N; a++){
		
		complex<double> z[N_TAU/2];
		
		for(int i=0; i<N_TAU/2; i++){
			double omega_n = (double)(2*i+1) * M_PI / beta;
// 			z[i] = IMAG * omega_n - G[a].self_c[i] + chem_pot;
			z[i] = IMAG * omega_n - G[a].self_c[i] + chem_pot - ec[a];
		}
		
		
		double ek_max=1.5*D;
		for(int ik=0; ik<N_p; ik++){
			ek[ik] = -ek_max + 2.0 * ek_max * (double)ik / (double)(N_p-1);
			
			ne[a][ik] = dmft_nk_cs_sub1(z, ek[ik], beta);
			ne_de[a][ik] = dmft_nk_cs_sub2(z, ek[ik], beta);
		}
		
		
		#ifdef DMFT_ISO
		mu_small[a] = eval_chem_pot(nc, D);
		mu_large[a] = eval_chem_pot(nc + 0.5, D);
		#else // DMFT_ISO
		mu_small[a] = eval_chem_pot(G[a].n_c, D);
		mu_large[a] = eval_chem_pot(G[a].n_c + 0.5, D);
		#endif // DMFT_ISO
		
		ne_small[a] = dmft_nk_cs_sub1(z, mu_small[a], beta);
		ne_large[a] = dmft_nk_cs_sub1(z, mu_large[a], beta);
		ne_de_small[a] = dmft_nk_cs_sub2(z, mu_small[a], beta);
		ne_de_large[a] = dmft_nk_cs_sub2(z, mu_large[a], beta);
		
		//
		// find local maxima of ne_de
		//  using the golden section search
		//
		n_peak[a]=0;
		int flag_pm=+1;
		
		for(int ik=1; ik<N_p; ik++){
			if( ne_de[a][ik+1] < ne_de[a][ik] ){
				
				if(flag_pm){
	// 				printf(" %.5e %.5e %.5e\n", ek[ik], ne[ik], ne_de[ik]);
					
					double x1=ek[ik-1], x2=ek[ik], x3=ek[ik+1];
					double y1=ne_de[a][ik-1], y2=ne_de[a][ik], y3=ne_de[a][ik];
					
					double golden = 0.5*(3.0-sqrt(5.0));
					
					do{
						if( x2-x1 > x3-x2 ){
							double x4 = x2 - golden * (x2-x1);
							double y4 = dmft_nk_cs_sub2(z, x4, beta);
							
							if(y4<y2){
								x1 = x4;
								y1 = y4;
							}
							else{
								x3 = x2;
								y3 = y2;
								x2 = x4;
								y2 = y4;
							}
							
						}
						else{
							double x4 = x2 + golden * (x3-x2);
							double y4 = dmft_nk_cs_sub2(z, x4, beta);
							
							if(y4<y2){
								x3 = x4;
								y3 = y4;
							}
							else{
								x1 = x2;
								y1 = y2;
								x2 = x4;
								y2 = y4;
							}
						}
						
					}while( x3-x1 > 1.0e-5 );
					
					ek_peak[a][n_peak[a]] = x2;
					ne_peak[a][n_peak[a]] = dmft_nk_cs_sub1(z, x2, beta);
					ne_de_peak[a][n_peak[a]] = y2;
					
					n_peak[a]++;
					if(n_peak[a]==3) break;
				}
				
				flag_pm = 0;
			}
			else{
				flag_pm = 1;
			}
		}
		
	}
	
	
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_ne.dat", num);
	
	FILE *fp=fopen(filename, "w");
	
	for(int a=0; a<N; a++){
		fprintf(fp, "# a=%d\n", a);
		fprintf(fp, "# small  %.6e %.6e %.6e\n", mu_small[a], ne_small[a], ne_de_small[a]);
		fprintf(fp, "# large  %.6e %.6e %.6e\n", mu_large[a], ne_large[a], ne_de_large[a]);
		
		for(int ik=0; ik<n_peak[a]; ik++){
			fprintf(fp, "# peak%d  %.6e %.6e %.6e\n", ik, ek_peak[a][ik], ne_peak[a][ik], ne_de_peak[a][ik]);
		}
	}
	
	for(int ik=0; ik<N_p; ik++){
		fprintf(fp, "%.6e", ek[ik]);
		for(int a=0; a<N; a++){
			fprintf(fp, " %.6e %.6e", ne[a][ik], ne_de[a][ik]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}


static void dmft_nk_ads_sub1(complex<double> *z, complex<double> *self_c, double V_sqr, double ek, double beta, double &nk_c, double &nk_f)
{
	// occupation number n(e)
	
	complex<double> Gk_c[N_TAU/2], Gk_f[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		Gk_c[i] = 1.0 / (z[i] - ek);
		Gk_f[i] = self_c[i] * (1.0 + self_c[i] * Gk_c[i]) / V_sqr;
	}
	
	double sum_Gk_c=0, sum_Gk_f=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk_c += real(Gk_c[i]);
		sum_Gk_f += real(Gk_f[i]);
	}
	sum_Gk_c *= 2.0 / beta;
	sum_Gk_f *= 2.0 / beta;
	
	sum_Gk_c += 0.5;
	sum_Gk_f += 0.5;
	
	// return
	nk_c = sum_Gk_c;
	nk_f = sum_Gk_f;
}
static void dmft_nk_ads_sub2(complex<double> *z, complex<double> *self_c, double V_sqr, double ek, double beta, double &nk_dk_c, double &nk_dk_f)
{
	// -dn(e)/de
	
	complex<double> Gk_c[N_TAU/2], temp[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		Gk_c[i] = 1.0 / (z[i] - ek);
		temp[i] = self_c[i] * Gk_c[i];
	}
	
	double sum_c=0, sum_f=0;
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_c += real(Gk_c[i] * Gk_c[i]) + 1.0 / (omega_n*omega_n);
		sum_f += real(temp[i] * temp[i])/V_sqr + 1.0 / (omega_n*omega_n);
	}
	sum_c *= 2.0 / beta;
	sum_f *= 2.0 / beta;
	
	sum_c -= 0.25 * beta;
	sum_f -= 0.25 * beta;
	
	// return
	nk_dk_c = -sum_c;
	nk_dk_f = -sum_f;
}
void dmft_nk_ads(complex<double> *self_c, double beta, double D, double V_sqr, double chem_pot, int num)
{
	
	complex<double> z[N_TAU/2];
	
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		z[i] = IMAG * omega_n - self_c[i] + chem_pot;
	}
	
	
	int N_p=1001;
	double ek[N_p], ne_c[N_p], ne_f[N_p], ne_de_c[N_p], ne_de_f[N_p];
	
	double ek_max=1.5*D;
	for(int ik=0; ik<N_p; ik++){
		ek[ik] = -ek_max + 2.0 * ek_max * (double)ik / (double)(N_p-1);
		
		dmft_nk_ads_sub1(z, self_c, V_sqr, ek[ik], beta, ne_c[ik], ne_f[ik]);
		dmft_nk_ads_sub2(z, self_c, V_sqr, ek[ik], beta, ne_de_c[ik], ne_de_f[ik]);
	}
	
	
// 	double mu_small = eval_chem_pot(nc, D);
// 	double mu_large = eval_chem_pot(nc + 0.5, D);
// 	
// 	double ne_small = dmft_nk_ads_sub1(z, mu_small, beta);
// 	double ne_large = dmft_nk_ads_sub1(z, mu_large, beta);
// 	double ne_de_small = dmft_nk_ads_sub2(z, mu_small, beta);
// 	double ne_de_large = dmft_nk_ads_sub2(z, mu_large, beta);
	
	
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_ne.dat", num);
	
	FILE *fp=fopen(filename, "w");
// 	fprintf(fp, "# small  %.6e %.6e %.6e\n", mu_small, ne_small, ne_de_small);
// 	fprintf(fp, "# large  %.6e %.6e %.6e\n", mu_large, ne_large, ne_de_large);
	
// 	for(int ik=0; ik<n_peak; ik++){
// 		fprintf(fp, "# peak%d  %.6e %.6e %.6e\n", ik, ek_peak[ik], ne_peak[ik], ne_de_peak[ik]);
// 	}
	for(int ik=0; ik<N_p; ik++){
		fprintf(fp, "%.6e %.6e %.6e %.6e %.6e\n", ek[ik], ne_c[ik], ne_de_c[ik], ne_f[ik], ne_de_f[ik]);
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
}


static void dmft_nk_ads_sub1(complex<double> *zf, complex<double> *zc, double V_sqr, double f_disp, double ek, double beta, double &nk_f, double &nk_c)
{
	// occupation number n(e)
	
	complex<double> Gk_c[N_TAU/2], Gk_f[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		Gk_f[i] = 1.0 / (zf[i] - f_disp*ek - V_sqr/(zc[i] - ek));
		Gk_c[i] = 1.0 / (zc[i] - ek - V_sqr/(zf[i] - f_disp*ek));
	}
	
	double sum_Gk_c=0, sum_Gk_f=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk_c += real(Gk_c[i]);
		sum_Gk_f += real(Gk_f[i]);
	}
	sum_Gk_c *= 2.0 / beta;
	sum_Gk_f *= 2.0 / beta;
	
	sum_Gk_c += 0.5;
	sum_Gk_f += 0.5;
	
	// return
	nk_c = sum_Gk_c;
	nk_f = sum_Gk_f;
}
static void dmft_nk_ads_sub2(complex<double> *zf, complex<double> *zc, double V_sqr, double f_disp, double ek, double beta, double &nk_dk_f, double &nk_dk_c)
{
	// -dn(e)/de
	
	nk_dk_c = 0;
	nk_dk_f = 0;
	
/*	complex<double> Gk_c[N_TAU/2], temp[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		Gk_c[i] = 1.0 / (z[i] - ek);
		temp[i] = self_c[i] * Gk_c[i];
	}
	
	double sum_c=0, sum_f=0;
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_c += real(Gk_c[i] * Gk_c[i]) + 1.0 / (omega_n*omega_n);
		sum_f += real(temp[i] * temp[i])/V_sqr + 1.0 / (omega_n*omega_n);
	}
	sum_c *= 2.0 / beta;
	sum_f *= 2.0 / beta;
	
	sum_c -= 0.25 * beta;
	sum_f -= 0.25 * beta;
	
	// return
	nk_dk_c = -sum_c;
	nk_dk_f = -sum_f;
	*/
	
}
void dmft_nk_ads(complex<double> **self_f, int N, double beta, double D, double *V_sqr, double *ef, double f_disp, double *ec, double chem_pot, int num)
{
	
	int N_p=1001;
	double ek[N_p], ne_c[N][N_p], ne_f[N][N_p], ne_de_c[N][N_p], ne_de_f[N][N_p];
	
	for(int s=0; s<N; s++){
		complex<double> zf[N_TAU/2], zc[N_TAU/2];
		
		for(int i=0; i<N_TAU/2; i++){
			double omega_n = (double)(2*i+1) * M_PI / beta;
	// 		z[i] = IMAG * omega_n - self_c[i] + chem_pot;
			zc[i] = IMAG * omega_n + chem_pot - ec[s];
			zf[i] = IMAG * omega_n + chem_pot - ef[s] - self_f[s][i];
		}
		
		
		double ek_max=1.5*D;
		for(int ik=0; ik<N_p; ik++){
			ek[ik] = -ek_max + 2.0 * ek_max * (double)ik / (double)(N_p-1);
			
			dmft_nk_ads_sub1(zf, zc, V_sqr[s], f_disp, ek[ik], beta, ne_f[s][ik], ne_c[s][ik]);
			dmft_nk_ads_sub2(zf, zc, V_sqr[s], f_disp, ek[ik], beta, ne_de_f[s][ik], ne_de_c[s][ik]);
		}
	}
	
	
// 	double mu_small = eval_chem_pot(nc, D);
// 	double mu_large = eval_chem_pot(nc + 0.5, D);
// 	
// 	double ne_small = dmft_nk_ads_sub1(z, mu_small, beta);
// 	double ne_large = dmft_nk_ads_sub1(z, mu_large, beta);
// 	double ne_de_small = dmft_nk_ads_sub2(z, mu_small, beta);
// 	double ne_de_large = dmft_nk_ads_sub2(z, mu_large, beta);
	
	
	char filename[40];
	sprintf(filename, DATA_DIR "%02d-dmft_ne.dat", num);
	
	FILE *fp=fopen(filename, "w");
// 	fprintf(fp, "# small  %.6e %.6e %.6e\n", mu_small, ne_small, ne_de_small);
// 	fprintf(fp, "# large  %.6e %.6e %.6e\n", mu_large, ne_large, ne_de_large);
	
// 	for(int ik=0; ik<n_peak; ik++){
// 		fprintf(fp, "# peak%d  %.6e %.6e %.6e\n", ik, ek_peak[ik], ne_peak[ik], ne_de_peak[ik]);
// 	}
	for(int ik=0; ik<N_p; ik++){
		fprintf(fp, "%.6e", ek[ik]);
		for(int s=0; s<N; s++){
			fprintf(fp, " %.6e %.6e %.6e %.6e", ne_c[s][ik], ne_de_c[s][ik], ne_f[s][ik], ne_de_f[s][ik]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
}


static double dmft_nk_hubbard_sub1(complex<double> *F, double ek, double beta, double jump)
{
	// occupation number n(e)
	
	complex<double> Gk[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Gk[i] = 1.0 / ( 1./ F[i] - ek);
	
	double sum_Gk=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk += real(Gk[i]);
	}
	sum_Gk *= 2.0 / beta;
	sum_Gk += 0.5 * jump;
	
	return( sum_Gk );
}
static double dmft_nk_hubbard_sub2(complex<double> *F, double ek, double beta, double jump)
{
	// -dn(e)/de
	
	complex<double> Gk[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Gk[i] = 1.0 / ( 1./ F[i] - ek);
	
	double sum_Gk_sqr=0;
	for(int i=0; i<N_TAU/2; i++){
		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk_sqr += real(Gk[i] * Gk[i]) + jump * jump / (omega_n*omega_n);
	}
	sum_Gk_sqr *= 2.0 / beta;
	sum_Gk_sqr -= 0.25 * beta * jump * jump;
	
	return( -sum_Gk_sqr );
}
static double dmft_nk_hubbard_sub3(complex<double> *F, double ek, double beta, double jump)
{
	// d^2 n(e)/de^2
	
	complex<double> Gk[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  Gk[i] = 1.0 / ( 1./ F[i] - ek);
	
	double sum_Gk=0;
	for(int i=0; i<N_TAU/2; i++){
// 		double omega_n = (double)(2*i+1) * M_PI / beta;
		sum_Gk += real(Gk[i] * Gk[i] * Gk[i]);
	}
	sum_Gk *= 2.0 / beta;
	
	return( sum_Gk );
}

struct params_dmft_nk_hubbard{
	complex<double> *F;
	double beta;
	double jump;
};
static double func_nk_hubbard(double ek, void *params)
{
	struct params_dmft_nk_hubbard *p = (params_dmft_nk_hubbard *)params;
	return dmft_nk_hubbard_sub3(p->F, ek, p->beta, p->jump);
}

// output: double ne_peak[s][3] = {ek, n(ek), -dn(ek)/ek}
void dmft_nk_hubbard(struct dmft_green_func *G, int N, double ne_peak[][3], double beta, double D, int num, double *jump)
{
	const int N_p=1001;
	double ek[N_p], ne[N][N_p], ne_de[N][N_p];
	double ne_de_2[N][N_p];
	double ek_max = D * 1.5;
	
	for(int s=0; s<N; s++){
		for(int ik=0; ik<N_p; ik++){
			ek[ik] = -ek_max + 2.0 * ek_max * (double)ik / (double)(N_p-1);
			
			ne[s][ik] = dmft_nk_hubbard_sub1(G[s].F, ek[ik], beta, jump[s]);
			ne_de[s][ik] = dmft_nk_hubbard_sub2(G[s].F, ek[ik], beta, jump[s]);
			ne_de_2[s][ik] = dmft_nk_hubbard_sub3(G[s].F, ek[ik], beta, jump[s]);
		}
	}
	
	//
	// find minimum of dn(e)/de by solving d^2 n(e)/de^2 = 0
	//
	for(int s=0; s<N; s++){
		double x_lo = -ek_max, x_hi = ek_max;
		struct params_dmft_nk_hubbard p = {G[s].F, beta, jump[s]};
		
		gsl_function F;
		F.function = &func_nk_hubbard;
		F.params = &p;
		
// 		printf(" x_lo = %.6lf\n", F.function(x_lo, F.params));
// 		printf(" x_hi = %.6lf\n", F.function(x_hi, F.params));
		
		if( F.function(x_lo, F.params) > 0 || F.function(x_hi, F.params) < 0 ){
			ne_peak[s][0] = ne_peak[s][1] = ne_peak[s][2] = NAN;
		} else {
	// 		gsl_root_fsolver *solver = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
			gsl_root_fsolver *solver = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
			gsl_root_fsolver_set (solver, &F, x_lo, x_hi);
			
	// 		printf ("using %s method\n", gsl_root_fsolver_name (solver));
	// 		printf ("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root", "err");
			
			int status;
			int iter = 0, max_iter = 100;
			double r = 0;
			do{
				iter++;
				status = gsl_root_fsolver_iterate (solver);
				r = gsl_root_fsolver_root (solver);
				x_lo = gsl_root_fsolver_x_lower (solver);
				x_hi = gsl_root_fsolver_x_upper (solver);
				status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7);
				
	// 			if (status == GSL_SUCCESS)  printf ("Converged:\n");
	// 			printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
			}
			while (status == GSL_CONTINUE && iter < max_iter);
			
			gsl_root_fsolver_free (solver);
			
			ne_peak[s][0] = r;
			ne_peak[s][1] = dmft_nk_hubbard_sub1(G[s].F, r, beta, jump[s]);
			ne_peak[s][2] = dmft_nk_hubbard_sub2(G[s].F, r, beta, jump[s]);
		}
	}
	
	
	char filename[128];
	sprintf(filename, DATA_DIR "%02d-dmft_ne.dat", num);
	FILE *fp=fopen(filename, "w");
	
	for(int a=0; a<N; a++){
		fprintf(fp, "# %.6e %.6e %.6e\n", ne_peak[a][0], ne_peak[a][1], ne_peak[a][2]);
	}
	for(int ik=0; ik<N_p; ik++){
		fprintf(fp, "%.6e", ek[ik]);
		for(int a=0; a<N; a++){
			fprintf(fp, " %.6e %.6e", ne[a][ik], ne_de[a][ik]);
			fprintf(fp, " %.6e", ne_de_2[a][ik]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}
void dmft_nk_hubbard(struct dmft_green_func *G, int N, double ne_peak[][3], double beta, double D, int num)
{
	double jump[N];
	for(int s=0; s<N; s++)  jump[s] = 1;
	dmft_nk_hubbard(G, N, ne_peak, beta, D, num, jump);
}


static double eval_nc(double mu)
{
	double nc=0;
	
	switch(DOS){
		case 0:
			if( fabs(mu) < 1.0 ){
				nc = 0.5 * (mu + 1.0);
			}
			else{
				if(mu > 0)  nc = 1;
				else        nc = 0;
			}
			break;
			
		case 1:
			if( fabs(mu) < 1.0 ){
				nc = 0.5 + ( mu * sqrt(1-mu*mu) + asin(mu) ) / M_PI;
			}
			else{
				if(mu > 0)  nc = 1;
				else        nc = 0;
			}
			break;
			
		case 2:
			nc = 0.5 * (1.0 + gsl_sf_erf(sqrt(2)*mu));
			break;
	}
	
	return(nc);
}

static double eval_chem_pot(double nc, double D)
{
	double chem;
	
	switch(DOS){
		case 0:
			chem = 2.0*nc - 1.0;
			break;
		case 1:
		case 2:
			double chem_max = 5;
			double d_chem = chem_max;
			double chem_accuracy = 1e-7;
			
			chem = chem_max;
			do{
				
				if( eval_nc(chem) - nc > 0 )  chem -= d_chem;
				else                      chem += d_chem;
				
				d_chem *= 0.5;
				
			}while( d_chem > chem_accuracy );
			
			break;
	}
	
	return(chem*D);
}



//=========================================================================================================
//
// thermodynamics
// (constant and symmetric DOS)
//
double internal_energy_imp(complex<double> *tmat, double beta, double D, double a)
{
	double E=0;
	
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		
		// subtract a/(i omega)
		double t_dif = imag(tmat[i]) + a / omega;
		
		E += 2.0/beta * omega / (omega*omega + D*D) * t_dif;
	}
	
	E += a / (2.0*D) * ( distrib_fermi(beta*D) - distrib_fermi(-beta*D) );
	
	return(E);
}


double internal_energy_perio(complex<double> *Gc, double beta, double D, double eff_mu)
{
	double E=0;
	
	// (3) subtract G0_omega
	complex<double> g0_omega[1][N_TAU/2];
	double temp_nc, E_c[1]={0};
	G0_omega_calc(g0_omega, 1, beta, D, 0, temp_nc, eff_mu, E_c);
	
	
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		
		// (1) subtract 1/(i omega)
// 		double Gc_dif = imag(Gc[i]) + 1.0 / omega;
		
		// (2) subtract 1/(i omega + eff_mu)
// 		double Gc_dif = imag( Gc[i] - 1.0 / (IMAG*omega + eff_mu) );
		
		// (3) subtract G0_omega
		double Gc_dif = imag(Gc[i]) - imag(g0_omega[0][i]);
		
		
		E -= omega * Gc_dif;
	}
	
	E *= 2.0 / beta;
	
	
	// (1)
// 	E -= 0.5 * eff_mu;
	
	// (2)
// 	double f = 0.5 * ( 1.0 - tanh(0.5*beta*(-eff_mu)) );
// 	E -= eff_mu * f;
	
	// (3)
	double eval_energy_free(double, double, double);
	E += eval_energy_free(beta, eff_mu, D);
	
	
	return(E);
}


double internal_energy_perio_ads(complex<double> *Gc, complex<double> *self_c, double beta, double D, double eff_mu)
{
	double E=0;
	
	// (3) subtract G0_omega
// 	complex<double> g0_omega[1][N_TAU/2];
// 	double temp_nc, E_c[1]={0};
// 	G0_omega_calc(g0_omega, 1, beta, D, 0, temp_nc, eff_mu, E_c);
	
	
	for(int i=0; i<N_TAU/2; i++){
		complex<double> i_omega = (double)(2*i+1) * M_PI / beta * IMAG;
		
		E += real( (i_omega + self_c[i]) * Gc[i] );
		
		// (1) subtract 1/(i omega)
// 		double Gc_dif = imag(Gc[i]) + 1.0 / omega;
		E -= 1.0;
		
		// (2) subtract 1/(i omega + eff_mu)
// 		double Gc_dif = imag( Gc[i] - 1.0 / (IMAG*omega + eff_mu) );
		
		// (3) subtract G0_omega
// 		double Gc_dif = imag(Gc[i]) - imag(g0_omega[0][i]);
// 		E -= real( i_omega * g0_omega[0][i] );
		
		
	}
	
	E *= 2.0 / beta;
	
	
	// (1)
	E -= 0.5 * eff_mu;
	
	// (2)
// 	double f = 0.5 * ( 1.0 - tanh(0.5*beta*(-eff_mu)) );
// 	E -= eff_mu * f;
	
	// (3)
// 	double eval_energy_free(double, double, double);
// 	E += eval_energy_free(beta, eff_mu, D);
	
	
	return(E);
}

// Evaluate E = T\sum_n iw G(iw)
// The coefficient a, b, c are required
//  a/iw, b/(iw)^2, c/(iw)^3  (normally a = 1, b = -mu)
double internal_energy_3(complex<double> *G, double beta, double a, double b, double c)
{
// 	printf("\n internal_energy_3()\n");
// 	printf(" a = %.4e\n", a);
// 	printf(" b = %.4e\n", b);
// 	printf(" c = %.4e\n", c);
// 	
// 	char filename[128];
// 	sprintf(filename, DATA_DIR "test-internal_energy.dat");
// 	FILE *fp=fopen(filename, "w");
// 	for(int i=0; i<N_TAU/2; i++){
// 		double omega = (double)(2*i+1) * M_PI / beta;
// 		fprintf(fp, "%.5e", omega);
// 		fprintf(fp, " %.6e", imag(G[i]));
// 		fprintf(fp, " %.6e", imag(G[i]) + a/omega);
// 		fprintf(fp, " %.6e", imag(G[i]) + a/omega - c/pow(omega,3));
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	
	double E=0;
	for(int i=0; i<N_TAU/2; i++){
		double omega = (double)(2*i+1) * M_PI / beta;
		
		// subtract a/(iw) and c/(iw)^3
		double G_dif = imag(G[i]) + a/omega - c/pow(omega,3);
		
		E -= omega * G_dif;
	}
	E *= 2. / beta;
	
	E += b / 2.;
	E -= c  * beta / 4.;
	
// 	printf(" E = %.5e\n", E);
	return(E);
}

static int sign(double x)
{
	return (x > 0) - (x < 0);
}
static int check_osc(int s1, int s2)
{
	if( s1 != 0 && s1 != s2 )  return 1;  // sign oscillates
	return 0;
}
// Estimate c and evaluate E
double internal_energy_2(complex<double> *G, double beta, double a, double b)
{
	// extrapolate y(x) to x=0
	// least-square method
// 	int num = 100;
	int num = MIN(100, N_TAU/2/10);
	double x2=0, x4=0, y1=0, y1x2=0;
	int sign_y = 0;
	bool flag_osc = 0;
	for(int i=0; i<num; i++){
		int j = N_TAU/2 - i - 1;
		double omega = (double)(2*j+1) * M_PI / beta;
		double x = 1. / omega;
		double y = ( imag(G[j]) + a/omega ) * pow(omega,3);
		
		y1 += y;
		y1x2 += y * pow(x,2);
		x2 += pow(x,2);
		x4 += pow(x,4);
		
		// check oscillation of y
		flag_osc += check_osc(sign_y, sign(y));
		sign_y = sign(y);
	}
	y1 /= double(num);
	y1x2 /= double(num);
	x2 /= double(num);
	x4 /= double(num);
	
	// average
	double c = y1;
	
	// quadratic function
// 	double c = ( y1 * x4 - y1x2 * x2 ) / (x4 - x2*x2 );
	
	// If oscillate, c does not make sense (the data is already converged within statistical errors)
	if( flag_osc )  c = 0;
	
	return internal_energy_3(G, beta, a, b, c);
}

double eval_energy_free(double beta, double eff_mu, double D)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	
	double result, error;
	double alpha[3] = {beta, eff_mu, D};
	
	double func_energy(double, void *);
	
	gsl_function F;
	F.function = func_energy;
	F.params = alpha;
	
	#if DOS==0 || DOS==1
	gsl_integration_qag(&F, -D, D, 0, 1e-7, 1000, 6, w, &result, &error);
	
	#elif DOS==2
	gsl_integration_qagi(&F, 0, 1e-7, 1000, w, &result, &error);
	
	#endif
	
	return(result);
}
double func_energy(double y, void *params)
{
	
	double *alpha = (double *)params;
	
	double beta = alpha[0];
	double eff_mu = alpha[1];
	double D = alpha[2];
	
	double dos_free(double, double);
	
	return( (y-eff_mu) * dos_free(y, D) * distrib_fermi(beta*(y-eff_mu)) );
	
}
double dos_free(double y, double D)
{
	
	double r = 0;
	
	#if DOS==0
	if( fabs(y) < D )  r = 0.5 / D;
	
	#elif DOS==1
	if( fabs(y) < D )  r = 2.0 / M_PI / D * sqrt( 1.0 - pow(y/D, 2) );
	
	#elif DOS==2
	r = sqrt( 2.0/M_PI ) / D * exp( -2.0 * pow(y/D, 2) );
	
	#endif
	
	return(r);
}


//=========================================================================================================




static complex<double> I_site_diag_2d(complex<double> pi, double J)
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			double Jq = J * 2.* ( cosk[i] + cosk[j] );
// 			double Jq = J[0] * 2.* ( cosk[i] + cosk[j] );
// 			 + J[1] * 4.* cosk[i] * cosk[j];
			
			r += w[i] * w[j] * Jq / (1. - Jq * pi);
		}
	}
	
	return( r / double(N_L*N_L) );
}

static complex<double> I_site_diag_3d(complex<double> pi, double J)
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				double Jq = J * 2.* ( cosk[i] + cosk[j] + cosk[k] );
// 				double Jq = J[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] )
// 				 + J[1] * 8.* cosk[i] * cosk[j] * cosk[k];
				
				r += w[i] * w[j] * w[k] * Jq / (1. - Jq * pi);
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L) );
}

static complex<double> I_site_diag_4d(complex<double> pi, double J)
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				for(int l=0; l<=N_L; l++){
					double Jq = J * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] );
// 					double Jq = J[0] * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] )
// 					 + J[1] * 16.* cosk[i] * cosk[j] * cosk[k] * cosk[l];
					
					r += w[i] * w[j] * w[k] * w[l] * Jq / (1. - Jq * pi);
				}
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L*N_L) );
}

static complex<double> I_site_diag_5d(complex<double> pi, double J)
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				for(int l=0; l<=N_L; l++){
					for(int m=0; m<=N_L; m++){
						double Jq = J * 2.* ( cosk[i] + cosk[j] + cosk[k] + cosk[l] + cosk[m] );
						
						r += w[i] * w[j] * w[k] * w[l] * w[m] * Jq / (1. - Jq * pi);
					}
				}
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L*N_L*N_L) );
}

static complex<double> I_site_diag_dummy(complex<double> pi, double J)
{
	printf("\n***error: I_site_diag_dummy()\n");
	exit(0);
}

// static complex<double> I_site_diag(complex<double> pi, double J[2])
// {
// 	#if DOS == 3
// 		return I_site_diag_2d(pi, J);
// 	#elif DOS == 4
// 		return I_site_diag_3d(pi, J);
// 	#elif DOS == 5
// 		return I_site_diag_4d(pi, J);
// 	#endif
// }

#if DOS == 3
	#define I_site_diag I_site_diag_2d
#elif DOS == 4
	#define I_site_diag I_site_diag_3d
#elif DOS == 5
	#define I_site_diag I_site_diag_4d
#elif DOS == 6
	#define I_site_diag I_site_diag_5d
#else
	#define I_site_diag I_site_diag_dummy  // dummy
#endif


// for test
static complex<double> chi_site_diag_3d(complex<double> pi1, complex<double> pi2, double J)
{
	double w[N_L+1];
	for(int i=1; i<N_L; i++)  w[i] = 1.0;
	w[0] = w[N_L] = 0.5;
	
	double cosk[N_L+1];
	for(int i=0; i<=N_L; i++){
		double k = double(i) * M_PI / double(N_L);
		cosk[i] = cos(k);
	}
	
	complex<double> temp = 1./ (pi1*pi2);
	complex<double> r = 0;
	for(int i=0; i<=N_L; i++){
		for(int j=0; j<=N_L; j++){
			for(int k=0; k<=N_L; k++){
				double Jq = J * 2.* ( cosk[i] + cosk[j] + cosk[k] );
				
				r += w[i] * w[j] * w[k] / (temp - Jq * Jq);
			}
		}
	}
	
	return( r / double(N_L*N_L*N_L) );
}

// J: +FM, -AFM
void sdmft1(sdmft_interaction &I, complex<double> *chi_sp, double J)
{// 
// 	double prm_J[2] = {J, 0};  // nearest neighbor
	
	for(int i=0; i<=N_TAU/2; i++){
		I.pi[i] = 1./ ( 1./ chi_sp[i] + I.I_eff[i] );
		
// 		I.chi_loc[i] = hilbert_transform_3d(1./ I.pi[i], prm_J);
// 		I.I_eff[i] = 1./ I.pi[i] - 1./ I.chi_loc[i];
		
		complex<double> I_diag = I_site_diag(I.pi[i], J);
		
		I.chi_loc[i] = I.pi[i] + I.pi[i] * I_diag * I.pi[i];
		
		I.I_eff[i] = 1./ ( 1./ I_diag + I.pi[i] );
	}
}



int sdmft_nearestneighbor()
{
	int z = 0;
	
	if( DOS >= 3 && DOS <= 6 ){
		int DIM = DOS - 1;
		// DIM = 2 for DOS = 3
		// DIM = 3 for DOS = 4
		// DIM = 4 for DOS = 5
		z = 2 * DIM;
	}
	else{
		printf("sdmft_nearestneighbor()\n");
		exit(0);
	}
	
	return z;
}

double sdmft_J_sqr(double J)
{
	double J_sqr = 0;
	
	if( DOS >= 3 && DOS <= 6 ){
		J_sqr = (double)sdmft_nearestneighbor() * J * J;
	}
	
	return J_sqr;
}


struct chi2pi{
	gsl_interp_accel *acc;
	gsl_spline *spline;
	double pi_max, chi_max;
};

static struct chi2pi *chi2pi_init()
{
	static chi2pi C2P;
	
	int N = 1001;
	
	C2P.spline = gsl_spline_alloc (gsl_interp_cspline, N);
	C2P.acc = gsl_interp_accel_alloc ();
	
// 	FILE *fp = NULL;
// 	if( DOS == 3 )  fp = fopen("chi_diag_2d.dat", "r");
// 	if( DOS == 4 )  fp = fopen("chi_diag_3d.dat", "r");
// 	if( fp == NULL ){
// 		printf("\n*** error chi2pi_init\n");
// 		exit(0);
// 	}
	FILE *fp;
	char filename[64];
	if( DOS == 3 )  sprintf(filename, "chi_diag_2d.dat");
	if( DOS == 4 )  sprintf(filename, "chi_diag_3d.dat");
	if( DOS == 5 )  sprintf(filename, "chi_diag_4d.dat");
	if( DOS == 6 )  sprintf(filename, "chi_diag_5d.dat");
	if( ( fp = fopen(filename, "r") ) == NULL ){
		printf("\n*** error : `%s' not found\n", filename);
		exit(0);
	}
	
	double pi[N], chi[N];
	for(int i=0; i<N; i++){
		fscanf(fp, "%lf %lf", &pi[i], &chi[i]);
	}
	fclose(fp);
	
	gsl_spline_init (C2P.spline, chi, pi, N);
	
	C2P.pi_max = pi[N-1];
	C2P.chi_max = chi[N-1];
	
	return &C2P;
}

static double chi2pi_eval_2d(double chi)
{
	static chi2pi *C2P = chi2pi_init();
	
	if( chi > C2P->chi_max ){
		double b = C2P->chi_max + log(1.0 - C2P->pi_max) / M_PI;
		return( 1.0 - exp( -(chi - b) * M_PI ) );
	}
	else{
		return gsl_spline_eval(C2P->spline, chi, C2P->acc);
	}
}

static double chi2pi_eval_3d(double chi)
{
	static chi2pi *C2P = chi2pi_init();
	
	if( chi > C2P->chi_max ){
		return 1;
// 		return NAN;
// 		exit(0);
	}
	else{
		return gsl_spline_eval(C2P->spline, chi, C2P->acc);
	}
}


#if DOS == 3
	#define chi2pi_eval chi2pi_eval_2d
#elif DOS == 4 || DOS == 5 || DOS == 6
	#define chi2pi_eval chi2pi_eval_3d
#else
	#define chi2pi_eval chi2pi_eval_3d  // dummy
#endif

// J: +FM, -AFM
// input: chi_sp (pure real), moment
void sdmft2(sdmft_interaction &I, complex<double> *chi_sp, double moment, double J)
{
	if( DOS == 1 ){  // Bethe lattice
		for(int i=0; i<=N_TAU/2; i++){
			I.pi[i] = chi_sp[i] / ( 1. + chi_sp[i] * I.I_eff[i] );
			I.chi_loc[i] = hilbert_transform(1./ I.pi[i], fabs(J));
			
			I.I_eff[i] = 0.25 * J * J * chi_sp[i];
		}
		I.mf = 0;
	}
	else if( DOS >= 3 && DOS <= 6 ){
		
		double z = sdmft_nearestneighbor();
		double zJ = z * fabs(J);
		
		for(int i=0; i<=N_TAU/2; i++){
			// dimensionless
			double chi_sp_zJ = real(chi_sp[i]) * zJ;
			if( chi_sp_zJ > 1e-4 ){
				// chi_sp -> pi
				I.pi[i] = chi2pi_eval( chi_sp_zJ );
				I.pi[i] /= zJ;  // restore dimension
				
				I.chi_loc[i] = 1./ (1./ I.pi[i] - I.I_eff[i]);
				
				// pi -> I_eff
				I.I_eff[i] = 1./ I.pi[i] - 1./ real(chi_sp[i]);
			}
			else{
				I.pi[i] = chi_sp[i];
				I.chi_loc[i] = chi_sp[i];
				I.I_eff[i] = I.pi[i] * sdmft_J_sqr(J);
			}
		}
		
		I.mf = moment * ( z * J - real(I.I_eff[0]) );
	}
	else{
		exit(0);
	}
}

// J: +FM, -AFM
// input: chi_sp (complex), moment
void sdmft3(sdmft_interaction &I, complex<double> *chi_sp, double moment, double J)
{
	for(int i=1; i<=N_TAU/2; i++){
		I.pi[i] = chi_sp[i] / ( 1. + chi_sp[i] * I.I_eff[i] );
		
// 		I.chi_loc[i] = hilbert_transform(1./ I.pi[i], J);
// 		I.I_eff[i] = 1./ I.pi[i] - 1./ I.chi_loc[i];
		
		complex<double> I_diag = I_site_diag(I.pi[i], J);
		I.chi_loc[i] = I.pi[i] + I.pi[i] * I_diag * I.pi[i];
		I.I_eff[i] = I_diag / ( 1. + I_diag * I.pi[i] );
	}
	
	double z = sdmft_nearestneighbor();
	
	{  // static part
		double zJ = z * fabs(J);
		
		for(int i=0; i<1; i++){
			// dimensionless
			double chi_sp_zJ = real(chi_sp[i]) * zJ;
			if( chi_sp_zJ > 1e-4 ){
				// chi_sp -> pi
				I.pi[i] = chi2pi_eval( chi_sp_zJ );
				I.pi[i] /= zJ;  // restore dimension
				
				I.chi_loc[i] = 1./ (1./ I.pi[i] - I.I_eff[i]);
				
				// pi -> I_eff
				I.I_eff[i] = 1./ I.pi[i] - 1./ real(chi_sp[i]);
			}
			else{
				I.pi[i] = chi_sp[i];
				I.chi_loc[i] = chi_sp[i];
				I.I_eff[i] = I.pi[i] * sdmft_J_sqr(J);
			}
		}
	}
	
	I.mf = moment * ( z * J - real(I.I_eff[0]) );
}

// J: +FM, -AFM
// input: chi_sp (complex), moment
void sdmft_twosubl(sdmft_interaction **I, complex<double> **chi_sp, double *moment, double J)
{
	for(int l=0; l<2; l++){
		for(int i=0; i<=N_TAU/2; i++){
			I[l]->pi[i] = chi_sp[l][i] / ( 1. + chi_sp[l][i] * I[l]->I_eff[i] );
		}
	}
	for(int i=0; i<=N_TAU/2; i++){
		complex<double> xi = sqrt( I[0]->pi[i] * I[1]->pi[i] );
		complex<double> I_diag = I_site_diag(xi, J);
		
		complex<double> temp1 = 1. + xi * I_diag;
		complex<double> temp2 = xi * I_diag / temp1;
		
		for(int l=0; l<2; l++){
			I[l]->chi_loc[i] = temp1 * I[l]->pi[i];
			I[l]->I_eff[i] = temp2 / I[l]->pi[i];
		}
	}
// 	for(int i=0; i<=N_TAU/2; i++){
// 		complex<double> chi_diag = chi_site_diag_3d(I[0]->pi[i], I[1]->pi[i], J);
// 		
// 		for(int l=0; l<2; l++){
// 			I[l]->chi_loc[i] = chi_diag / I[(l+1)%2]->pi[i];
// 		}
// 		for(int l=0; l<2; l++){
// 			I[l]->I_eff[i] = 1./ I[l]->pi[i] - 1./ I[l]->chi_loc[i];
// 		}
// 	}
	
	double z = sdmft_nearestneighbor();
	
	{  // static part
		double zJ = z * fabs(J);
		
		for(int i=0; i<1; i++){
			// dimensionless
			double chi_sp_zJ = sqrt( real(chi_sp[0][i]) * real(chi_sp[1][i]) ) * zJ;
			if( chi_sp_zJ > 1e-4 ){
				// chi_sp -> pi
				double xi = chi2pi_eval( chi_sp_zJ );
				double r = sqrt( real(chi_sp[0][i]) / real(chi_sp[1][i]) );
				I[0]->pi[i] = xi * r;
				I[1]->pi[i] = xi / r;
				
				for(int l=0; l<2; l++){
					I[l]->pi[i] /= zJ;  // restore dimension
					
					I[l]->chi_loc[i] = 1./ (1./ I[l]->pi[i] - I[l]->I_eff[i]);
					
					// pi -> I_eff
					I[l]->I_eff[i] = 1./ I[l]->pi[i] - 1./ real(chi_sp[l][i]);
				}
			}
			else{
				for(int l=0; l<2; l++){
					I[l]->pi[i] = chi_sp[l][i];
					I[l]->chi_loc[i] = chi_sp[l][i];
					I[l]->I_eff[i] = I[l]->pi[i] * sdmft_J_sqr(J);
				}
			}
		}
	}
	
	for(int l=0; l<2; l++){
		I[l]->mf = moment[(l+1)%2] * z * J - moment[l] * real(I[l]->I_eff[0]);
	}
}



static void print_sdmft_sub1(sdmft_interaction *I, int N, char *filename)
{
	FILE *fp=fopen(filename, "w");
	
	for(int i=0; i<=N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int a=0; a<N; a++){
			fprintf(fp, " %.8e", real(I[a].I_eff[i]));
			fprintf(fp, " %.8e", real(I[a].chi_loc[i]));
			fprintf(fp, " %.8e", real(I[a].pi[i]));
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	printf("\n '%s'\n", filename);
}

void print_sdmft(sdmft_interaction *I, int N, int num, int subl)
{
	char filename[40];
	if( subl < 0 ){
		sprintf(filename, DATA_DIR "%02d-sdmft.dat", num);
	}
	else{
		sprintf(filename, DATA_DIR "%02d-%c-sdmft.dat", num, subl+'A');
	}
	
	print_sdmft_sub1(I, N, filename);
}

void print_sdmft(sdmft_interaction *I, int N, int num)
{
	print_sdmft(I, N, num, -1);
}


static void print_sdmft_complex_sub1(sdmft_interaction *I, int N, char *filename)
{
	FILE *fp=fopen(filename, "w");
	
	for(int i=0; i<=N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int a=0; a<N; a++){
			fprintf(fp, " %.8e %.8e", real(I[a].I_eff[i]), imag(I[a].I_eff[i]));
			fprintf(fp, " %.8e %.8e", real(I[a].chi_loc[i]), imag(I[a].chi_loc[i]));
			fprintf(fp, " %.8e %.8e", real(I[a].pi[i]), imag(I[a].pi[i]));
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	printf("\n '%s'\n", filename);
}

void print_sdmft_complex(sdmft_interaction *I, int N, int num, int subl)
{
	char filename[40];
	if( subl < 0 ){
		sprintf(filename, DATA_DIR "%02d-sdmft.dat", num);
	}
	else{
		sprintf(filename, DATA_DIR "%02d-%c-sdmft.dat", num, subl+'A');
	}
	
	print_sdmft_complex_sub1(I, N, filename);
}

void print_sdmft_complex(sdmft_interaction *I, int N, int num)
{
	print_sdmft_complex(I, N, num, -1);
}


// 1: I is set
// 0: I is not set
int read_sdmft(struct sdmft_interaction *I, int N, int num, int subl)
{
	FILE *fp;
	char filename[64], prefix[64];
	
	if(subl<0)  sprintf(prefix, "%02d", num);
	else        sprintf(prefix, "%02d-%c", num, subl+'A');
	
	sprintf(filename, "%s-sdmft.dat", prefix);
	
	if( (fp=fopen(filename, "r")) != NULL ){
		
// 		if(my_rank==0){
// 			printf("\n---\n");
// 			printf("\n '%s' opened\n", filename);
// 		}
		
		int n=0;  // data number
		double data[128];
		
		while( read_data(fp, data) >= 1+3*N ){
			for(int s=0; s<N; s++){
				I[s].I_eff[n] = data[1+3*s];
				I[s].chi_loc[n] = data[2+3*s];
				I[s].pi[n] = data[3+3*s];
			}
			n++;
		}
		fclose(fp);
		
		// test
// 		for(int i=0; i<n; i++){
// 			printf("%d", i);
// 			for(int a=0; a<N_F; a++)  printf(" %.5e %.5e", real(G0_omega[a][i]), imag(G0_omega[a][i]));
// 			printf("\n");
// 		}
		
		
		if(n != N_TAU/2 + 1){
// 			if(my_rank==0){
				printf("\n*** error : inconsistency in data number\n");
				printf("  n = %d,  N_TAU/2 = %d\n", n, N_TAU/2);
// 			}
			exit(0);
		}
// 		if(my_rank==0)  printf(" G0_omega[] has been set\n");
		return 1;
	}
	else{
// 		if(my_rank==0){
// 			printf("\n---\n");
// 			printf("\n '%s' does not exist\n", filename);
// 		}
		return 0;
	}
}

int read_sdmft(struct sdmft_interaction *I, int N, int num)
{
	return read_sdmft(I, N, num, -1);
}


// 1: I is set
// 0: I is not set
int read_sdmft_complex(struct sdmft_interaction *I, int N, int num, int subl)
{
	FILE *fp;
	char filename[64], prefix[64];
	
	if(subl<0)  sprintf(prefix, "%02d", num);
	else        sprintf(prefix, "%02d-%c", num, subl+'A');
	
	sprintf(filename, "%s-sdmft.dat", prefix);
	
	if( (fp=fopen(filename, "r")) != NULL ){
		
// 		if(my_rank==0){
// 			printf("\n---\n");
// 			printf("\n '%s' opened\n", filename);
// 		}
		
		int n=0;  // data number
		double data[128];
		
		while( read_data(fp, data) >= 1+6*N ){
			for(int s=0; s<N; s++){
				I[s].I_eff[n]   = complex<double>(data[1+6*s], data[2+6*s]);
				I[s].chi_loc[n] = complex<double>(data[3+6*s], data[4+6*s]);
				I[s].pi[n]      = complex<double>(data[5+6*s], data[6+6*s]);
			}
			n++;
		}
		fclose(fp);
		
		// test
// 		for(int i=0; i<n; i++){
// 			printf("%d", i);
// 			for(int a=0; a<N_F; a++)  printf(" %.5e %.5e", real(G0_omega[a][i]), imag(G0_omega[a][i]));
// 			printf("\n");
// 		}
		
		
		if(n != N_TAU/2 + 1){
// 			if(my_rank==0){
				printf("\n*** error : inconsistency in data number\n");
				printf("  n = %d,  N_TAU/2 = %d\n", n, N_TAU/2);
// 			}
			exit(0);
		}
// 		if(my_rank==0)  printf(" G0_omega[] has been set\n");
		return 1;
	}
	else{
// 		if(my_rank==0){
// 			printf("\n---\n");
// 			printf("\n '%s' does not exist\n", filename);
// 		}
		return 0;
	}
}

int read_sdmft_complex(struct sdmft_interaction *I, int N, int num)
{
	return read_sdmft_complex(I, N, num, -1);
}
