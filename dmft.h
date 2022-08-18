/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _DMFT
#define _DMFT

#include "ct_qmc_share.h"
// dependences:
// N_TAU
// DOS
// ACCURACY_FILLING
// FILLING_ITER_MAX
// DATA_DIR
//
#include "cpplapack_suppl.h"


//
// Green functions for the periodic systems (DMFT)
//
struct dmft_green_func{
// 	double Gc_cav_tau[N_TAU];
	complex<double> Gc_cav_omega[N_TAU/2];
	complex<double> self_c[N_TAU/2], *F;
	complex<double> Gc_perio[N_TAU/2];
	double n_c;
	complex<double> tmat_perio[N_TAU/2], *Gf_perio;
	double n_f;  // only for Anderson lattice
// 	complex<double> chi_c_unif[N_TAU/2], chi_c_stag[N_TAU/2];
// 	complex<double> T0_unif[N_TAU/2], T0_stag[N_TAU/2];
	dmft_green_func(){
		F = self_c;
		Gf_perio = tmat_perio;
	}
};

struct dmft_tp{
	complex<double> Pi_loc[N_TAU/2], Pi_unif[N_TAU/2], Pi_stag[N_TAU/2];
	complex<double> T0_loc[N_TAU/2], T0_unif[N_TAU/2], T0_stag[N_TAU/2];
	complex<double> T0_loc2[N_TAU/2];
	double stat_Pi_loc, stat_Pi_unif, stat_Pi_stag;
};


// Hilbert transform
// output:
//     (1/N) sum_k (z - ek)^{-1} = int dos(e) / (z - e)
// input: z
//        D (half bandwidth or n.n. hopping)
// 
complex<double> hilbert_transform(complex<double> z, double D);

// evaluate G(tau=-0), occupation number
double eval_G_tau0m(complex<double> *G_omega, double beta);
double eval_G_tau0m(complex<double> *G_omega, double beta, double jump);


// DMFT for CS and Kondo lattice
void dmft_cs(struct dmft_green_func &G, complex<double> *tmat, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot);
void dmft_cs(struct dmft_green_func *G, complex<double> **tmat, int N, double beta, double D, int flag_ensemble, double &n_c, double &chem_pot, double *E_c);

void dmft_cs_2sub(struct dmft_green_func *G[2], complex<double> **tmat[2], int N, double beta, double D, int flag_ensemble, double &ave_n_c, double &chem_pot, double *E_c);

// DMFT for Anderson lattice

// f-level fix
// input: Gf
void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr);
// n_f, tot_n: particle numbers per spin
void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot);
void dmft_ads0(struct dmft_green_func &G, complex<double> *Gf, double beta, double D, double V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double jump);
// n_f, tot_n: particle numbers (not per spin)
void dmft_ads0(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double *E_c, double *jump);

void dmft_ads0_2sub(struct dmft_green_func *G[2], complex<double> **Gf[2], int N, double beta, double D, double *V_sqr, double n_f, int flag_ensemble, double &tot_n, double &chem_pot, double *E_c, double *jump[2]);

double dmft_ads_Veff(double f_disp, double D);
double dmft_Veff(double D);

// input: self_f
void dmft_ads1(struct dmft_green_func *G, complex<double> **self_f, int N, double beta, double D, double *V_sqr, double *ef, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c);

// input: Gf
void dmft_ads2(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c, double *jump);
void dmft_ads2(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, double *V_sqr, double f_disp, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c);

// input: tmat, nf, Gf
// for CPA with hybridization disorder
void dmft_ads3(struct dmft_green_func *G, complex<double> **tmat, double *nf_in, int N, double beta, double D, int flag_ensemble, double &ave_n, double &chem_pot, double *E_c);

// input: Gf
void dmft_hubbard(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, int flag_ensemble, double &occup, double &chem_pot, double *jump);
void dmft_hubbard(struct dmft_green_func *G, complex<double> **Gf, int N, double beta, double D, int flag_ensemble, double &occup, double &chem_pot);

// layer Hubbard model
// input: Gf
void dmft_layer_hb(struct dmft_green_func **G, complex<double> **Gf, int N, double beta, double D, const CPPL::dsymatrix &hop);
void dmft_layer_current(struct dmft_green_func **G, int N, double beta, double D, const CPPL::dsymatrix &hop, double *Kz, double **Kxy, int N_nu);


void dmft_two_particle(struct dmft_green_func &G, struct dmft_tp &T0, double beta, double D, int n_f, double chem_pot);
void dmft_two_particle_dynamic(struct dmft_green_func &G, complex<double> *Pi_unif_dynamic, double beta, double D, double chem_pot);

void print_dmft_num(struct dmft_green_func *G, int N, int num1, int num2);
void print_dmft_num(struct dmft_green_func *G, int N, double pot_scat, int num1, int num2);

void print_dmft(struct dmft_green_func *G, int N, double pot_scat, int num, int subl);
void print_dmft(struct dmft_green_func *G, int N, double pot_scat, int num);
void print_dmft(struct dmft_green_func *G, int N, int num, int subl);
void print_dmft(struct dmft_green_func *G, int N, int num);

void print_dmft_mu(struct dmft_green_func *G, int N, int num, int subl, double mu);
void print_dmft_mu(struct dmft_green_func *G, int N, int num, double mu);

int read_dmft(struct dmft_green_func *G, int N, int num, int subl);
int read_dmft(struct dmft_green_func *G, int N, int num);

int read_dmft_mu(struct dmft_green_func *G, int N, int num, int subl, double &mu);
int read_dmft_mu(struct dmft_green_func *G, int N, int num, double &mu);

void print_dmft_tp(struct dmft_tp &T0, int num);
void print_dmft_tp_dynamic(complex<double> *Pi_unif_dynamic, double beta, double D, int num);


void dmft_nk_cs(struct dmft_green_func &G, double beta, double D, double chem_pot, double nc, int num);
void dmft_nk_cs(struct dmft_green_func *G, int N, double beta, double D, double chem_pot, double nc, int num);
void dmft_nk_cs(struct dmft_green_func *G, int N, double beta, double D, double chem_pot, double nc, double *ec, int num);

void dmft_nk_ads(complex<double> *self_c, double beta, double D, double V_sqr, double chem_pot, int num);
void dmft_nk_ads(complex<double> **self_f, int N, double beta, double D, double *V_sqr, double *ef, double f_disp, double *ec, double chem_pot, int num);

// output: double ne_peak[s][3] = {ek, n(ek), -dn(ek)/ek}
void dmft_nk_hubbard(struct dmft_green_func *G, int N, double ne_peak[][3], double beta, double D, int num, double *jump);
void dmft_nk_hubbard(struct dmft_green_func *G, int N, double ne_peak[][3], double beta, double D, int num);



//
// thermodynamics
//
double internal_energy_imp(complex<double> *tmat, double beta, double D, double a);
double internal_energy_perio(complex<double> *Gc, double beta, double D, double eff_mu);
double internal_energy_perio_ads(complex<double> *Gc, complex<double> *self_c, double beta, double D, double eff_mu);

// Evaluate E = T\sum_n iw G(iw)
// The coefficient a, b, c are required
//  a/iw, b/(iw)^2, c/(iw)^3  (normally a = 1, b = -mu)
double internal_energy_3(complex<double> *G, double beta, double a, double b, double c);

// Estimate c and evaluate E
double internal_energy_2(complex<double> *G, double beta, double a, double b);



//
// DMFT for exchange interaction
//

struct sdmft_interaction{
	double mf;  // molecular field
	complex<double> I_eff[N_TAU/2+1];  // +: FM
	complex<double> chi_loc[N_TAU/2+1];
	complex<double> pi[N_TAU/2+1];
};


// J: +FM, -AFM
void sdmft1(sdmft_interaction &I, complex<double> *chi_sp, double J);
// chi_sp (real)
void sdmft2(sdmft_interaction &I, complex<double> *chi_sp, double moment, double J);
// chi_sp (complex)
void sdmft3(sdmft_interaction &I, complex<double> *chi_sp, double moment, double J);
// two-sublattice
void sdmft_twosubl(sdmft_interaction **I, complex<double> **chi_sp, double *moment, double J);

int sdmft_nearestneighbor();
double sdmft_J_sqr(double J);

void print_sdmft(sdmft_interaction *I, int N, int num, int subl);
void print_sdmft(sdmft_interaction *I, int N, int num);
void print_sdmft_complex(sdmft_interaction *I, int N, int num, int subl);
void print_sdmft_complex(sdmft_interaction *I, int N, int num);

int read_sdmft(struct sdmft_interaction *I, int N, int num, int subl);
int read_sdmft(struct sdmft_interaction *I, int N, int num);
int read_sdmft_complex(struct sdmft_interaction *I, int N, int num, int subl);
int read_sdmft_complex(struct sdmft_interaction *I, int N, int num);


#endif // _DMFT
