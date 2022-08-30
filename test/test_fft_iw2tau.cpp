#include "fft.h"
#include <cstdio>
#include <cmath>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

main()
{
    double beta=10;
    int N=1024;
    double e0 = 0.4;
    double tail = 0.77;
    double acc = 1e-3;

    vector< complex<double> > G_iw(N/2);
    for(int i=0; i<N/2; i++){
        double omega = double(2*i+1) * M_PI / beta;
        G_iw[i] = tail / (IMAG*omega - e0);
    }

    vector<double> G_tau;
    fft_fermion_iw2tau(G_tau, G_iw, beta, tail);

    double error = 0;
    printf("# size of G_tau: %d\n", G_tau.size());
    for(int i=0; i<N; i++){
        double tau = beta * double(i) / double(N);
        double G_ref = -tail * exp(-e0*tau) / ( exp(-e0*beta) + 1.0 );
        error = max( error, fabs(G_ref - G_tau[i]) );
        printf("%d %lf %lf\n", i, G_tau[i], G_ref);
    }

    printf(" max error = %.2e\n", error);
    if( error > acc ){
        printf(" ERROR\n");
    }
}
