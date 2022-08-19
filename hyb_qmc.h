/*

Continuous-time quantum Monte Carlo
for the impurity Anderson model (hybridization expansion)

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _HYB_QMC_H
#define _HYB_QMC_H

#include "ct_qmc_share.h"
#include "operators.h"
#include <vector>
#include "vector_type.hpp"
#include "array.hpp"

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
// #define CHI_TR 0 // transverse susceptibility (N_S==2)

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

// struct single_particle{
// 	double Gf_tau[N_TAU+1], Gf_tau_err[N_TAU+1];
// 	complex<double> Gf_omega[N_TAU/2];
// 	double f_number, f_number_err;
// 	double f_number2, f_number2_err;
// 	double jump;  // coefficient of 1/iw tail; jump=1 for UINF=0, jump<1 for UINF=1
// // 	complex<double> Gf0_omega[N_TAU/2];
// 	complex<double> self_f[N_TAU/2];
// };
struct single_particle{
	vec_d Gf_tau, Gf_tau_err;
	vec_c Gf_omega;
	double f_number, f_number_err;
	double f_number2, f_number2_err;
	double jump;  // coefficient of 1/iw tail; jump=1 for UINF=0, jump<1 for UINF=1
// 	complex<double> Gf0_omega[N_TAU/2];
	vec_c self_f;

	single_particle() {};
	single_particle(int n_tau){
		Gf_tau.resize(n_tau+1);
		Gf_tau_err.resize(n_tau+1);
		Gf_omega.resize(n_tau/2);
		self_f.resize(n_tau/2);
	};
};

// struct two_particle{
// 	double chi_tau[N_TP+1], chi_tau_err[N_TP+1];
// 	complex<double> chi_omega[N_TP2+1];
// };
struct two_particle{
	vec_d chi_tau, chi_tau_err;
	vec_c chi_omega;

	two_particle() {};
	two_particle(int n_tp, int n_tp2){
		chi_tau.resize(n_tp+1);
		chi_tau_err.resize(n_tp+1);
		chi_omega.resize(n_tp2+1);
	};
};

// transverse susceptibility
// struct two_particle_tr{
// 	double chi_tau1[2*N_TP2+1], chi_tau1_err[2*N_TP2+1];
// 	double chi_tau2[2*N_TP2+1], chi_tau2_err[2*N_TP2+1];
// 	complex<double> chi_omega1[N_TP2+1], chi_omega2[N_TP2+1];
// 	complex<double> chi_omega[N_TP_W], chi_omega_err[N_TP_W];  // w-sampling
// };

// struct phys_quant{
// 	double ave_sign, ave_sign_err;
// 	double Z_k[N_S][N_K], Z_k_err[N_S][N_K];
// 	double Z_ktot[N_S*N_K], Z_ktot_err[N_S*N_K];
// 	double ave_k[N_S], ave_k_err[N_S];
// 	double ave_ktot, ave_ktot_err;
// 	double occup_tot, occup_tot_err;
// 	double occup_mom, occup_mom_err;
// 	double stat_suscep_sp, stat_suscep_sp_err;
// 	double stat_suscep_ch, stat_suscep_ch_err;
// };
struct phys_quant{
	double ave_sign, ave_sign_err;
	vec_vec_d Z_k, Z_k_err;
	vec_d Z_ktot, Z_ktot_err;
	vec_d ave_k, ave_k_err;
	double ave_ktot, ave_ktot_err;
	double occup_tot, occup_tot_err;
	double occup_mom, occup_mom_err;
	double stat_suscep_sp, stat_suscep_sp_err;
	double stat_suscep_ch, stat_suscep_ch_err;

	phys_quant() {};
	phys_quant(int n_s, int n_k){
		resize(Z_k, n_s, n_k);
		resize(Z_k_err, n_s, n_k);
		Z_ktot.resize(n_s*n_k);
		Z_ktot_err.resize(n_s*n_k);
	};
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

//============================================================================

class HybQMC{
public:
	// initializing : allocating memory and creating log file
	HybQMC(int max_k, int n_s, int n_tau, int n_tp, int n_tp2, int rand_seed);
	~HybQMC();

	void get(struct phys_quant **p_PQ, struct single_particle **p_SP, double **p_TP_tau, struct two_particle ***p_TP, struct two_particle **p_TP_sp, struct two_particle **p_TP_ch, struct phonon_g **D);

	// setting number of MC steps
	void set_nmc(struct num_mc n_mc_in);

	// setting parameter set
	void set_params(struct hyb_qmc_params prm_in);

	// setting hybridization function
	//  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
	// void set_Delta(complex<double> Delta_omega[N_S][N_TAU/2], double V_sqr[N_S]);
	void set_Delta(vec_vec_c &Delta_omega, vec_d &V_sqr);
	//  Delta = V_sqr * G0
	// void set_G0(complex<double> G0_omega[N_S][N_TAU/2], double V_sqr[N_S]);

	// [Optional] setting moment, which is used for susceptibility calculations
	void set_moment(double moment_f_in[N_S]);

	// evaluating physical quantities via MC sampling
	//  flag_tp: 0 only single-particle, 1 single- and two-particle
	void eval(int flag_tp);

	// free memory
	// void final();

	// printing a string in log file
	void fprint_log(char *str);

private:

	struct hyb_qmc_params prm;
	struct num_mc n_mc;

	// static struct cond_op S[N_S];
	std::vector<Operators> S;
	std::vector<struct green_func_0> Delta;

	std::vector<struct single_particle> SP;
	std::vector<std::vector<struct two_particle> > TP;
	struct two_particle TP_sp, TP_ch;
	// struct two_particle_tr TP_tr;
	struct phys_quant PQ;
	struct phonon_g D;

	// static double TP_tau[N_TP+1];
	vec_d TP_tau;

	// static complex<double> Delta_omega[N_S][N_TAU/2];
	vec_vec_c Delta_omega;

	// static double moment_f[N_S];
	vec_d moment_f;

	static int w_sign;

	void eval_acceptance(double n_sample, int n_add, int n_shift);
	void sampling(int i_measure, int n_bin, int n_sample, int n_add, int n_shift);
	void opt_n_mc(double accept_seg, double accept_shift);

	void init_mc_config();
	void init_measure();
	void init_measure_bin();

	void measure_stat();
	void measure_sp();
	void measure_tp();
	void func_measure0();
	void func_measure1();
	void func_measure2();
	// void (*func_measure[3])() = {func_measure0, func_measure1, func_measure2};

	void averagebin_stat(int n_sample);
	void averagebin_sp(int n_sample);
	void averagebin_tp(int n_sample);
	void func_averagebin0(int);
	void func_averagebin1(int);
	void func_averagebin2(int);
	// void (*func_averagebin[3])(int) = {func_averagebin0, func_averagebin1, func_averagebin2};

	void fft_chi_after_interp(vec_d &chi_tau, vec_c &chi_omega);
	void average_stat(int n_bin);
	void average_sp(int n_bin);
	void average_tp(int n_bin);
	void func_average0(int);
	void func_average1(int);
	void func_average2(int);
	// void (*func_average[3])(int) = {func_average0, func_average1, func_average2};

	double mpi_reduce_bin(int);
	double mpi_reduce_accept();


	void state0_change(int sigma);
	void state0_change_2(int sigma1, int sigma2);

	void add_seg(int sigma, int anti);
	void rem_seg(int sigma, int anti);
	// void (* func_add_remove[2])(int, int)={
	// 	add_seg, rem_seg
	// };

	void shift_tau1(int sigma, int i_tau1);
	void shift_tau2(int sigma, int i_tau2);
	// void (* func_shift_tau[2])(int, int)={
	// 	shift_tau1, shift_tau2
	// };

	struct phys_quant_bin{
		long int ave_sign;
		Array2D<unsigned long> n_k;  // [N_S][N_K]
		std::vector<unsigned long> n_ktot;  // [N_S*N_K]
		Array2D<double> Gf;  // [N_S][N_TAU+1]
		std::vector<double> f_number;  // [N_S]
		std::vector<int> f_number_int;  // [N_S]
		double occup_tot, occup_mom;
		Array3D<double> chi; // [N_S][N_S][N_TP+1]

		phys_quant_bin() {};
		phys_quant_bin(int n_k, int n_s, int n_tau, int n_tp) : n_k(n_s, n_k), n_ktot(n_s*n_k), Gf(n_s, n_tau+1), f_number(n_s), f_number_int(n_s), chi(n_s, n_s, n_tp+1) {};

		void allzeros();
	};
	phys_quant_bin B, B_TOT;
};

//============================================================================
//
// FUNCTIONS TO BE CALLED BY MAIN FUNCTION
//

// // initializing : allocating memory and creating log file
// void hybqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, double **p_TP_tau, struct two_particle ***p_TP, struct two_particle **p_TP_sp, struct two_particle **p_TP_ch, struct two_particle_tr **TP_tr, struct phonon_g **D);

// // setting number of MC steps
// void hybqmc_set_nmc(struct num_mc n_mc_in);

// // setting parameter set
// void hybqmc_set_params(struct hyb_qmc_params prm_in);

// // setting hybridization function
// //  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
// void hybqmc_set_Delta(complex<double> Delta_omega[N_S][N_TAU/2], double V_sqr[N_S]);
// //  Delta = V_sqr * G0
// void hybqmc_set_G0(complex<double> G0_omega[N_S][N_TAU/2], double V_sqr[N_S]);

// // [Optional] setting moment, which is used for susceptibility calculations
// void hybqmc_set_moment(double moment_f_in[N_S]);

// // evaluating physical quantities via MC sampling
// //  flag_tp: 0 only single-particle, 1 single- and two-particle
// void hybqmc_eval(int flag_tp);

// // free memory
// void hybqmc_final();

// // printing a string in log file
// void hybqmc_fprint_log(char *str);


#endif // _HYB_QMC_H
