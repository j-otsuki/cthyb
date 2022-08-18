/*

Pade approximation

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "pade.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;


void pade_init(complex<double> *z, complex<double> *u, complex<double> *a, int N)
{
// 	printf("init_pade...\n");
	
// 	complex<double> g[N][N];  // g_n(z_i)
	
	complex<double>* g[N];  // g_n(z_i)
	for(int i=0; i<N; i++){
		g[i] = (complex<double> *)malloc(N*sizeof(complex<double>));
		if(g[i]==NULL){
			printf("*** memory not allocated\n");
			exit(0);
		}
	}
	
	for(int i=0; i<N; i++) g[0][i] = u[i];
	for(int n=1; n<N; n++){
		for(int i=n; i<N; i++)  g[n][i] = ( g[n-1][n-1] / g[n-1][i] - 1.0 ) / (z[i] - z[n-1]);
	}
	
	for(int i=0; i<N; i++)  a[i] = g[i][i];
	
	for(int i=0; i<N; i++)  free(g[i]);
}

complex<double> pade_calc(complex<double> *z, complex<double> *a, complex<double> w, int N)
{
	complex<double> A0 = 0.0, A1 = a[0], A2, B0 = 1.0, B1 = 1.0, B2;
	
	for(int i=1; i<N; i++){  // (N-1)/2
		A2 = A1 + ( w - z[i-1] ) * a[i] * A0;
		B2 = B1 + ( w - z[i-1] ) * a[i] * B0;
		
// 		printf("%d  %e  %e\n", i, abs(A2), abs(B2));
		A1/=B2;
		A2/=B2;
		B1/=B2;
		B2/=B2;
		
		A0 = A1;
		A1 = A2;
		B0 = B1;
		B1 = B2;
	}
	
	return( A2 / B2 );
}
