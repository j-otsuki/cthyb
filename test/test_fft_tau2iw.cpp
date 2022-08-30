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
    double acc = 1e-5;

    vector<double> G_tau(N);
    for(int i=0; i<N; i++){
        double tau = beta * double(i) / double(N);
        G_tau[i] = -tail * exp(-e0*tau) / ( exp(-e0*beta) + 1.0 );
    }

    vector< complex<double> > G_iw;
    fft_fermion_tau2iw(G_tau, G_iw, beta, tail);

    double error = 0;
    printf("# size of G_iw: %d\n", G_iw.size());
    for(int i=0; i<N/2; i++){
        double omega = double(2*i+1) * M_PI / beta;
        complex<double> G_ref = tail / (IMAG*omega - e0);
        error = max( error, abs(G_ref - G_iw[i]) );
        printf("%d %lf %lf %lf %lf\n", i, G_iw[i].real(), G_iw[i].imag(), G_ref.real(), G_ref.imag());
    }

    printf(" max error = %.2e\n", error);
    if( error > acc ){
        printf(" ERROR\n");
    }
}
