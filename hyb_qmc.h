/*

Continuous-time quantum Monte Carlo
for the impurity Anderson model (hybridization expansion)

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _HYB_QMC_H
#define _HYB_QMC_H

#include "ct_qmc_share.h"
//
// N_TAU
// N_K
// struct cond_op
//

#define HYB_QMC_MPI 1
#define PHONON 0  // Holstein model

#define N_S 2  // number of local states
#define N_TP 32  // two-particle, non-linear mesh used in measurement of two-particle responce function
#define N_TP2 256  // two-particle, linear mesh to perform FFT
#define N_TP_W 10  // two-particle, w-sampling

#define TEST_MODE 0
#define DISPLAY 1
#define CHI_TR 0 // transverse susceptibility (N_S==2)

#define LOG_FILE DATA_DIR "hyb_qmc.log"

const unsigned long RAND_SEED = 0;  // 0: random (time),  >0: fixed number
const int N_WARMUP = 1000000;
const int MAX_R_CORR = 0;  // >1: maximum of corr_fac for opt_n_mc,  0: no correction
const int K_TOT_MIN = 20 * N_S;  // >1: minimum of k_tot used in opt_n_mc,  0: no correction


#if HYB_QMC_MPI
#include <mpi.h>
#endif // HYB_QMC_MPI


/* INPUT
struct hyb_qmc_params prm;
struct num_mc n_mc;
complex<double> G0_omega[N_TAU/2];
*/

/* OUTPUT
struct phys_quant *PQ;
struct single_particle *SP;  // [N_S]
struct two_particle **TP;  // [N_S][N_S]
struct two_particle *TP_sp, *TP_ch;
struct two_particle_tr *TP_tr;
double *TP_tau;  // [N_TP+1]
*/

//============================================================================
//
// FUNCTIONS TO BE CALLED BY MAIN FUNCTION
//

// initializing : allocating memory and creating log file
void hybqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, double **p_TP_tau, struct two_particle ***p_TP, struct two_particle **p_TP_sp, struct two_particle **p_TP_ch, struct two_particle_tr **TP_tr, struct phonon_g **D);

// setting number of MC steps
void hybqmc_set_nmc(struct num_mc n_mc_in);

// setting parameter set
void hybqmc_set_params(struct hyb_qmc_params prm_in);

// setting hybridization function
//  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
void hybqmc_set_Delta(complex<double> Delta_omega[N_S][N_TAU/2], double V_sqr[N_S]);
//  Delta = V_sqr * G0
void hybqmc_set_G0(complex<double> G0_omega[N_S][N_TAU/2], double V_sqr[N_S]);

// [Optional] setting moment, which is used for susceptibility calculations
void hybqmc_set_moment(double moment_f_in[N_S]);

// evaluating physical quantities via MC sampling
//  flag_tp: 0 only single-particle, 1 single- and two-particle
void hybqmc_eval(int flag_tp);

// free memory
void hybqmc_final();

// printing a string in log file
void hybqmc_fprint_log(char *str);


//============================================================================
//
// DEFINITION OF STRUCT
//

struct hyb_qmc_params{
// 	double V_sqr;
	double ef[N_S];
	double U[N_S][N_S];  // symmetric, diagonal is not used
	double beta;
	int UINF;
	double w_ph;  // phonon frequency
	double g;  // electron-phonon coupling
	
	hyb_qmc_params(){
		UINF = 0;
	}
};

struct single_particle{
	double Gf_tau[N_TAU+1], Gf_tau_err[N_TAU+1];
	complex<double> Gf_omega[N_TAU/2];
	double f_number, f_number_err;
	double f_number2, f_number2_err;
	double jump;  // coefficient of 1/iw tail; jump=1 for UINF=0, jump<1 for UINF=1
// 	complex<double> Gf0_omega[N_TAU/2];
	complex<double> self_f[N_TAU/2];
};

struct two_particle{
	double chi_tau[N_TP+1], chi_tau_err[N_TP+1];
	complex<double> chi_omega[N_TP2+1];
};

// transverse susceptibility
struct two_particle_tr{
	double chi_tau1[2*N_TP2+1], chi_tau1_err[2*N_TP2+1];
	double chi_tau2[2*N_TP2+1], chi_tau2_err[2*N_TP2+1];
	complex<double> chi_omega1[N_TP2+1], chi_omega2[N_TP2+1];
	complex<double> chi_omega[N_TP_W], chi_omega_err[N_TP_W];  // w-sampling
};

struct phys_quant{
	double ave_sign, ave_sign_err;
	double Z_k[N_S][N_K], Z_k_err[N_S][N_K];
	double Z_ktot[N_S*N_K], Z_ktot_err[N_S*N_K];
	double ave_k[N_S], ave_k_err[N_S];
	double ave_ktot, ave_ktot_err;
	double occup_tot, occup_tot_err;
	double occup_mom, occup_mom_err;
	double stat_suscep_sp, stat_suscep_sp_err;
	double stat_suscep_ch, stat_suscep_ch_err;
};

struct phonon_g{
	complex<double> d_omega[N_TP2];
	double occup;
};

struct num_mc{
// 	int N_WARMUP;
	int N_MSR;
	int N_BIN;
	int N_ADD;  // to be optimized if N_ADD<0.  R_ADD = -N_ADD / 10.  
	int N_SHIFT;  // to be optimized if N_SHIFT<0.  R_SHIFT = -N_SHIFT / 10.
// The most part of the configuration is expected to be updated,
// when R_ADD=1 (N_ADD=-10) or R_SHIFT=1 (N_SHIFT=-10). 
};

#endif // _HYB_QMC_H
