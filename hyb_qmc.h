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

#define TEST_MODE 0
#define DISPLAY 1

#define LOG_FILE "hyb_qmc.log"

const int N_WARMUP = 1000000;

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
// DEFINITION OF STRUCT (INPUT)
//

struct hyb_qmc_params{
	vec_d ef;  // [N_S]
	vec_vec_d U;  // [N_S][N_S], symmetric, diagonal is not used
	double beta;
	int UINF;

	hyb_qmc_params(){
		UINF = 0;
	}
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
//
// DEFINITION OF STRUCT (OUTPUT)
//

struct single_particle{
	vec_d Gf_tau, Gf_tau_err;  // [N_TAU+1]
	vec_c Gf_omega;  // [N_TAU/2];
	double f_number, f_number_err;
	double f_number2, f_number2_err;
	double jump;  // coefficient of 1/iw tail; jump=1 for UINF=0, jump<1 for UINF=1
	vec_c self_f;  // [N_TAU/2]

	single_particle() {};
	single_particle(int n_tau);
	void allzeros();
};

struct two_particle{
	vec_d chi_tau, chi_tau_err;  // [N_TP+1]
	vec_c chi_omega;  // [N_TP2+1]

	two_particle() {};
	two_particle(int n_tp, int n_tp2);
	void allzeros();
};

struct phys_quant{
	double ave_sign, ave_sign_err;
	vec_vec_d Z_k, Z_k_err;  // [N_S][N_K]
	vec_d Z_ktot, Z_ktot_err;  // [N_S*N_K]
	vec_d ave_k, ave_k_err;  // [N_S]
	double ave_ktot, ave_ktot_err;
	double occup_tot, occup_tot_err;
	double occup_mom, occup_mom_err;
	double stat_suscep_sp, stat_suscep_sp_err;
	double stat_suscep_ch, stat_suscep_ch_err;

	phys_quant() {};
	phys_quant(int n_s, int n_k);
	void allzeros();
};

//============================================================================

using t_sp = std::vector<struct single_particle>;
using t_tp = std::vector<std::vector<struct two_particle> >;

class HybQMC{
public:
	// initializing : allocating memory and creating log file
	HybQMC(int max_order, int n_s, int n_tau, int n_tp, int n_tp2, int rand_seed);
	~HybQMC();

	// setting number of MC steps
	void set_nmc(const num_mc& n_mc_in);

	// setting parameter set
	void set_params(const hyb_qmc_params& prm_in);

	// setting hybridization function
	//  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
	void set_Delta(const vec_vec_c& Delta_omega, const vec_d& V_sqr);
	//  Delta = V_sqr * G0

	// [Optional] setting moment, which is used for susceptibility calculations
	void set_moment(const vec_d& moment_f_in);

	// evaluating physical quantities via MC sampling
	//  flag_tp = false : only single-particle
	//            true  : single- and two-particle
	void eval(bool flag_tp);

	// get
	t_sp& get_SP() { return SP; };
	t_tp& get_TP() { return TP; };
	vec_d& get_TP_tau() { return TP_tau; };
	two_particle& get_TP_sp() { return TP_sp; };
	two_particle& get_TP_ch() { return TP_ch; };
	phys_quant& get_PQ() { return PQ; };

	// printing a string in log file
	void fprint_log(char *str);

private:
	int N_TAU;
	int N_K;
	int N_S;
	int N_TP, N_TP2;

	struct hyb_qmc_params prm;
	struct num_mc n_mc;

	// static struct cond_op S[N_S];
	std::vector<Operators> S;
	std::vector<struct green_func_0> Delta;

	t_sp SP;
	t_tp TP;
	vec_d TP_tau;
	two_particle TP_sp, TP_ch;
	phys_quant PQ;
	vec_vec_c Delta_omega;
	vec_d moment_f;

	int w_sign;

	// If N_ADD<0 and/or N_SHIFT<0,
	// N_ADD and/or N_SHIFT is optimized according to
	//   N_ADD = (ave_tot_k / acceptance) * R_ADD * corr_fac
	//   N_SHIFT = (2 * ave_tot_k / acceptance) * R_SHIFT * corr_fac
	// The most part of the configuration is expected to be updated, when R_ADD=1 or R_SHIFT=1.
	double R_ADD;  // set by R_ADD = -N_ADD / 10
	double R_SHIFT;  // set by R_SHIFT = -N_SHIFT / 10
	int N_ADD_MIN, N_SHIFT_MIN;  // N_S, 1
	// correction factor: [1:MAX_R_CORR]  determined by P(k)

	int MAX_R_CORR;  // >1: maximum of corr_fac for opt_n_mc,  0: no correction
	int K_TOT_MIN;  // >1: minimum of k_tot used in opt_n_mc,  0: no correction

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

	void averagebin_stat(int n_sample);
	void averagebin_sp(int n_sample);
	void averagebin_tp(int n_sample);
	void func_averagebin0(int);
	void func_averagebin1(int);
	void func_averagebin2(int);

	void fft_chi_after_interp(const vec_d& chi_tau, vec_c& chi_omega);
	void average_stat(int n_bin);
	void average_sp(int n_bin);
	void average_tp(int n_bin);
	void func_average0(int);
	void func_average1(int);
	void func_average2(int);

	double mpi_reduce_bin(int);
	double mpi_reduce_accept();

	void state0_change(int sigma);
	void state0_change_2(int sigma1, int sigma2);

	void add_seg(int sigma, int anti);
	void rem_seg(int sigma, int anti);

	void shift_tau1(int sigma, int i_tau1);
	void shift_tau2(int sigma, int i_tau2);

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
		phys_quant_bin(int n_k, int n_s, int n_tau, int n_tp);
		void allzeros();
	};
	phys_quant_bin B, B_TOT;
};

#endif // _HYB_QMC_H
