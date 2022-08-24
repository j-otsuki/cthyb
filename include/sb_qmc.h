/*

Continuous-time quantum Monte Carlo
for the Anderson model with spin-boson coupling

written by Junya Otsuki, 2012
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _SB_QMC_H
#define _SB_QMC_H

#include "ct_qmc_share.h"
// dependences:
// N_TAU
// N_K
//
// #include "matrix_update.h"
// #include "green_func_0.h"

#define QMC_MPI 1
// #define LIMIT_PARAM 0
// 0: g != 0,  V != 0
// 1: g_xy = 0  ( Anderson model + g_z )
// 2: V = 0  ( spin-boson model )
#define UINF 0
// 0: U finite
// 1: U infinity

#define N_TP 64  // two-particle, non-linear mesh used in measurement of two-particle responce function
#define N_TP2 N_TAU/2  // two-particle, linear mesh to perform FFT

#define N_VX1 10  // number of fermionic frequency of vertex
#define N_VX2 10  // number of bosonic frequency of vertex

#define VX_NFFT 1  // How to compute vertex
// 0: no FFT (direct sum)
// 1: Non-equidistant FFT
// 2: TEST: compare the above two
const int M_OS=4;  // oversampling factor

#define IMPROV_ESTIM 1  // use improved estimators for vertex calculations

#define SELF_VERTEX_ISO 1
// In computing self-energy and vertex with improved estimators, assume...
//  0: g_xy=0
//  1: g_xy=g_z (multiply factor 3)

// const double MIN_AVE_SIGN = 0.1;

#define TEST_MODE 0  // 1: error check,  2: error check & output config
#define DISPLAY 1

#define LOG_FILE DATA_DIR "sb_qmc.log"

const unsigned long RAND_SEED = 0;  // 0: random (time),  >0: fixed number
const int N_WARMUP = 1000000;
const int MAX_AVE_PERT_V = N_K;
const int MAX_AVE_PERT_G = N_K;



#if QMC_MPI
#include <mpi.h>
#endif // QMC_MPI



//============================================================================
//
// FUNCTIONS TO BE CALLED BY MAIN FUNCTION
//

void sbqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, struct two_particle **p_TP, struct vertex **p_VX);
void sbqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, struct two_particle **p_TP);
// void sbqmc_init();

void sbqmc_set_nmc(struct num_mc n_mc_in);

// 0: successful exit
// 1: ave_sign is lower than MIN_AVE_SIGN
int sbqmc_eval(int flag_tp);

// evaluate free energy by Wang-Landau sampling
int sbqmc_eval_thermo();

void sbqmc_final();

void sbqmc_fprint_log(char *str);

// Delta_omega_in[2][N_TAU/2]
//  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
// D0_omega[3][N_TAU/2+1]   {ch, sp(zz), sp(+-)}
//  fac_int[a]:
//   S_int = (1/2)   fac_int[a] \int sigma^0(tau) D0_omega[0](tau-tau') sigma^0(tau')
//         + (1/2)   fac_int[a] \int sigma^z(tau) D0_omega[1](tau-tau') sigma^z(tau')
//         + (1/2) 4 fac_int[a] \int sigma^-(tau) D0_omega[2](tau-tau') sigma^+(tau')
// fac_chi[a]:
//  chi[0] =   fac_chi[0] < sigma^0 ; sigma^0 >
//  chi[1] =   fac_chi[1] < sigma^z ; sigma^z >
//  chi[2] = 2 fac_chi[2] < sigma^+ ; sigma^- >
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double *fac_int, double *fac_chi_in, double D_jump);
// D_jump = 0
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double *fac_int, double *fac_chi_in);
// D_jump = 0
// fac_int[3] = {1,1,1};
// fac_chi[3] = {1,1,1};
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in);


void sbqmc_set_input1(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double D_jump);
// D_jump = 0
void sbqmc_set_input1(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in);

// interactions are defined with the Pauli operator:  g * sigma . phi
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double D_jump);
// D_jump = 0
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in);

// // Delta = V_sqr * G0
// void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> G0_omega[][N_TAU/2], double V_sqr[], complex<double> **D0_omega_in, double D_jump);
// // D_jump = 0
// void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> G0_omega[][N_TAU/2], double V_sqr[], complex<double> **D0_omega_in);

void sbqmc_set_output(struct statistics *p_STAT, struct phys_quant *p_PQ, struct single_particle *p_SP, struct two_particle *p_TP, struct vertex *p_VX);


//============================================================================
//
// DEFINITION OF STRUCT
//

struct sb_qmc_params{
	double beta;
	double ef[2];
	double U;
};

// struct statistics{
// 	double d_accept_seg, d_reject_seg;
// 	double d_accept_spin, d_reject_spin;
// 	double d_accept_boson, d_reject_boson;
// 	double d_accept_shift, d_reject_shift;
// // 	double n_k[2][N_K];  // V
// // 	double n_l[N_K];  // g
// // 	unsigned long tot_n_k[N_A*N_K];
// };

struct single_particle{
	double Gf_tau[N_TAU+1], Gf_tau_err[N_TAU+1];
	complex<double> Gf_omega[N_TAU/2];
	double jump;
	complex<double> self_f_1[N_TAU/2];  // ordinary definition (not converge)
	complex<double> self_f_2[N_TAU/2];  // converge ~w^{-1}
	double G_self[N_TAU+1];  // G*Sigma
	complex<double> self_f[N_TAU/2];
};

struct two_particle{
	double tau[N_TP+1];
	double chi_sp_tau[N_TP+1], chi_sp_tau_err[N_TP+1];
	double chi_ch_tau[N_TP+1], chi_ch_tau_err[N_TP+1];
	complex<double> chi_sp_omega[N_TP2+1];
	complex<double> chi_ch_omega[N_TP2+1];
	
	double chi_pm_tau[N_TAU+1], chi_pm_tau_err[N_TAU+1];
	complex<double> chi_pm_omega[N_TAU/2+1];
};

struct vertex{
	complex<double> G[2][2*N_VX1+N_VX2];
	complex<double> self[2][2*N_VX1+N_VX2];  // set SELF_VERTEX_ISO
	complex<double> suscep_lo[2][2][N_VX2];
	// (omega, omega'; nu>=0)
	// *_lo : longitudinal [s1][s2]
	// *_tr : transverse   [2] = {+-, -+}
// 	complex<double> chi_sp[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> chi_ch[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> chi_pm[2*N_VX1][2*N_VX1][N_VX2];
	complex<double> Gfour_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];  // Two-particle Green function
	complex<double> Gfour_tr[2][2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx2_sp[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx2_ch[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx2_pm[2*N_VX1][2*N_VX1][N_VX2];
	complex<double> gamma2_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];  // four-point vertex
	complex<double> gamma2_tr[2][2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx_sp[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx_ch[2*N_VX1][2*N_VX1][N_VX2];
// 	complex<double> vx_pm[2*N_VX1][2*N_VX1][N_VX2];
	complex<double> gamma_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];  // set SELF_VERTEX_ISO
	complex<double> gamma_tr[2][2*N_VX1][2*N_VX1][N_VX2];  // set SELF_VERTEX_ISO
	// (omega; nu>=0)
	complex<double> Gthree_lo[2][2][2*N_VX1][N_VX2];
	complex<double> lambda_lo[2][2][2*N_VX1][N_VX2];  // three-point vertex
	complex<double> lambda_lo_test[2][2][2*N_VX1][N_VX2];
};

// struct spin_boson{
// 	complex<double> D_z_omega[N_TAU/2+1];
// 	complex<double> D_xy_omega[N_TAU/2+1];
// };

struct phys_quant{
	double ave_sign, ave_sign_err;
	double Z_V[2][N_K], Z_V_err[2][N_K];  // V
	double Z_g[N_K], Z_g_err[N_K];  // g
	double ave_Z_V[2], ave_Z_V_err[2];
	double ave_Z_g, ave_Z_g_err;
	double occup_tot, occup_tot_err;
	double occup_dif, occup_dif_err;
	double occup_ud[2], occup_ud_err[2];
	double occup_dbl, occup_dbl_err;
	double stat_suscep_sp, stat_suscep_sp_err;
	double stat_suscep_ch, stat_suscep_ch_err;
	double stat_suscep_pm, stat_suscep_pm_err;
	double m_z_sqr, m_z_sqr_err;
	double m_z_pw4, m_z_pw4_err;
	double m_xy_sqr, m_xy_sqr_err;
	double m_xy_pw4, m_xy_pw4_err;
	double m_rt_sqr, m_rt_sqr_err;
	double m_rt_pw4, m_rt_pw4_err;
	double eff_mom, stat_suscep_reg;
	double free_energy, free_energy_err;
};

struct num_mc{
// 	int WARMUP;
	int msr;
	int bin;
	int seg;
	int spin;
	int boson;
	int shift;
	int segsp;
};


#endif // _SB_QMC_H
