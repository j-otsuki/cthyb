/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _GREEN_FUNC_0_H
#define _GREEN_FUNC_0_H

#include "ct_qmc_share.h"
// dependences:
// N_TAU
// DOS
// ACCURACY_FILLING
// FILLING_ITER_MAX


double G0_calc_integ(double tau, double beta, double D);

void G0_omega_calc(std::complex<double> G0_omega[][N_TAU/2], int N, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c);
// void G0_omega_calc(std::complex<double> *G0_omega, double beta, double D, int flag_dmft, double &n_c, double &chem_pot);
void G0_omega_calc(std::complex<double> *G0_omega, double beta, double D);
void G0_omega_calc(std::complex<double> *G0_omega, double beta, double D, double V_sqr);

void Gf_calc_fft(double *Gf_tau, double beta, double D, double epsilon_f, double U, double V_sqr);


// DOS = delta(w - w0)
void D0_omega_calc0(std::complex<double> *D0_omega, double beta, double w0);

// DOS \propto omega^{gamma}  (with cutoff)
void D0_omega_calc1(std::complex<double> *D0_omega, double beta, double cutoff, double gamma);


#endif // _GREEN_FUNC_0_H
