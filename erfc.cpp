#include "erfc.h"
#include <stdio.h>
#include <math.h>
using namespace std;

static const double ERFC_ACCURACY=1.0e-12;

static inline complex<double> erfc_w_low(complex<double> z)
{
	complex<double> iz(-imag(z), real(z));
	complex<double> izn = 1;
	complex<double> w=1;
	
	double fac_odd = 1./sqrt(M_PI);  // n=1
	double fac_even = 1;  // n=2
	
	int i=1;
	complex<double> dw;
	
	do{
		fac_odd *= 2./(double)(2*i-1);
		izn *= iz;
		dw = fac_odd * izn;  // n=2*i-1
		
		fac_even /= (double)i;
		izn *= iz;
		dw += fac_even * izn;  // n=2*i
		
		w += dw;
		i++;
		
	}while( abs(dw) > ERFC_ACCURACY );
	
// 	printf("# i=%d\n", i);
	
	complex<double> fac(0, -sqrt(M_PI));
	
// 	return(w);  // w(z)
	return(fac * w);
}

static inline complex<double> erfc_w_high(complex<double> z)
{
// 	complex<double> A0=1, A1=0, A2;
// 	complex<double> B0=0, B1=1, B2;
	
	complex<double> A0=0, A1=1, A2;
	complex<double> B0=1, B1=z, B2;
	
	int i=1;
	double abs_A0, abs_A1=1;
	
	do{
		double a = -0.5 * (double)i;
		
		A2 = z * A1 + a * A0;
		B2 = z * B1 + a * B0;
		
		A0 = A1;
		B0 = B1;
		A1 = A2;
		B1 = B2;
		
		A0 /= B1;
		A1 /= B1;
		B0 /= B1;
		B1 = 1.0;
		
		abs_A0 = abs_A1;
		abs_A1 = abs(A1);
		
// 		printf(" %d", i);
		i++;
		
		if(i>1000)  break;
		
	}while( fabs(abs_A1 - abs_A0) > ERFC_ACCURACY );
	
// 	printf("# i=%d\n", i);
	
// 	complex<double> fac(0, 1./sqrt(M_PI));
	
// 	return(fac * A1);  // w(z)
	return(A1);
}

complex<double> erfc_w(complex<double> z)
{
	complex<double> w;
	
	if( imag(z) < 2.0 - 0.5 * real(z) )  w = erfc_w_low(z);
	else  w = erfc_w_high(z);
	
	return(w);
}

/*
main()  // for test
{
	
// 	complex<double> z0(1, 0);
// 	complex<double> z0(0, 1);
	complex<double> z0(cos(M_PI/8.), sin(M_PI/8.));
// 	complex<double> z0(1, 0.01);
	
	for(int i=0; i<70; i++){
		double r = 0. + 7. * (double)i / 70.;
		
		
// 		complex<double> w1 = erfc_w_low(z0 * r);
// 		complex<double> w2 = erfc_w_high(z0 * r);
// 		
// 		printf("%.3e %.8e %.8e %.8e %.8e\n", r, real(w1), imag(w1), real(w2), imag(w2));
		
		
		complex<double> w = erfc_w(z0 * r);
		
		printf("%.3e %.8e %.8e\n", r, real(w), imag(w));
	}
	
}
*/
