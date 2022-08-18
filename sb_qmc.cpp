/*

Continuous-time quantum Monte Carlo
for the Anderson model with spin-boson coupling

written by Junya Otsuki, 2012
Dept. of Physics, Tohoku University, Sendai, Japan

*/

static const char tag[64] = "v1.10";

#include "sb_qmc.h"
#include "matrix_update.h"
// #include "green_func_0.h"
#include <math.h>
#include "float.h"  // for DBL_MAX
#include "pade.h"
#include "common.h"
#include <fftw3.h>

static struct sb_qmc_params prm;
static struct num_mc n_mc;

static struct cond_op S[2];
static struct green_func_0 Delta[2];
// static struct boson_func_0 D0;  // g^2 * D0
// static struct boson_func_0 K0, K1;
// static struct boson_func_0 K0_ch, K1_ch;
// static struct boson_func_0 K0_sp, K1_sp;
static struct boson_func_0 K0[2];  // ch, sp
static struct boson_func_0 K1[2];
static struct boson_func_0 D0_pm;

static struct phys_quant PQ;
static struct single_particle *SP;
static struct two_particle *TP;
static struct vertex *VX;

static complex<double> Delta_omega[2][N_TAU/2];
// static complex<double> D0_xy_omega[N_TAU/2+1],  D0_z_omega[N_TAU/2+1];
static complex<double> D0_omega[3][N_TAU/2+1];  // g^2 * D0  [ch, sp, pm]
static int flag_V;
static int flag_g[3];
static double fac_chi[3];

static double prm_ef_unrenorm[2];
static double prm_U_unrenorm;
// static double TP_tau[N_TP+1];

static int w_sign;
static struct x_op{
	int k;
	double tau[4*N_K+1];
	int op[4*N_K+1];  // 0-5 (fu-, fu+, fd-, fd+, S-, S+)
	int state[4*N_K+1];  // state before the operator is acted
} X;
static struct b_op{
	int k;
	double tau1[N_K], tau2[N_K];
} Y;

enum OP { F_UP_AN, F_UP_CR, F_DW_AN, F_DW_CR, S_MINUS, S_PLUS};
enum OP F_OP[2][2] = {  // [sigma][dag]
	{F_UP_AN, F_UP_CR},
	{F_DW_AN, F_DW_CR}
};
enum OP S_OP[2] = {S_MINUS, S_PLUS};

// enum FSTATE {UP, DW, EMP};
// enum FSTATE op2state_after[6] = {EMP, UP, EMP, DW, DW, UP}; // The state AFTER the operator acts
// enum FSTATE op2state[6] = {UP, EMP, DW, EMP, UP, DW}; // The state BEFORE the operator acts
enum FSTATE {EMP, UP, DW, DBL};  // 00, 01, 10, 11

// static double s_fac_ch[6] = {-1.0,  1.0, -1.0,  1.0,  0., 0.};
// static double s_fac_sp[6] = {-0.5,  0.5,  0.5, -0.5, -1., 1.};
static double s_fac[2][6] = {
	{-1.0,  1.0, -1.0,  1.0,  0., 0.},
	{-0.5,  0.5,  0.5, -0.5, -1., 1.}
};
// static double *s_fac_ch = s_fac[0];
// static double *s_fac_sp = s_fac[1];

struct mc_accept{
	double d_accept_seg, d_reject_seg;
	double d_accept_spin, d_reject_spin;
	double d_accept_boson, d_reject_boson;
	double d_accept_shift, d_reject_shift;
	double d_accept_segsp, d_reject_segsp;
// 	double n_k[2][N_K];  // V
// 	double n_l[N_K];  // g
// 	unsigned long tot_n_k[N_A*N_K];
} ACCPT;
static struct mc_accept_sub{
	unsigned long n_accept_seg, n_reject_seg;
	unsigned long n_accept_spin, n_reject_spin;
	unsigned long n_accept_boson, n_reject_boson;
	unsigned long n_accept_shift, n_reject_shift;
	unsigned long n_accept_segsp, n_reject_segsp;
// 	unsigned long n_k[2][N_K];  // V
// 	unsigned long n_l[N_K];  // g
} ACCPT_sub;

static struct phys_quant_bin{
	long int ave_sign;
	unsigned long n_k[2][N_K];  // V
	unsigned long n_l[N_K];  // g
	double Gf[2][N_TAU+1];
	double self[2][N_TAU+1];
	double occup[4];
	double m_z_sqr, m_z_pw4;
	double m_xy_sqr, m_xy_pw4;
	double m_rt_sqr, m_rt_pw4;
	double chi_pm[N_TAU+1];
	double chi[4][4][N_TP+1];
// 	double chi_tr1[2*N_TP2+1], chi_tr2[2*N_TP2+1];
	complex<double> vx_Gf[2][2*N_VX1+N_VX2];
	complex<double> vx_self[2][2*N_VX1+N_VX2];
	complex<double> vx_G4_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];
	complex<double> vx_gamma_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];
	complex<double> vx_G4_tr[2][2*N_VX1][2*N_VX1][N_VX2];
	complex<double> vx_gamma_tr[2][2*N_VX1][2*N_VX1][N_VX2];
	complex<double> vx_suscep[2][2][N_VX2];
	complex<double> vx_G3[2][2][2*N_VX1][N_VX2];
} *B, *B_TOT;

static struct vertex_temp{
	complex<double> phase_i[2*N_VX1+N_VX2][N_K];
	complex<double> phase_j[2*N_VX1+N_VX2][N_K];
	complex<double> u[2*N_VX1+N_VX2][2*N_VX1+N_VX2];
	complex<double> u2[2*N_VX1+N_VX2][2*N_VX1+N_VX2];
	int flag_occup1[N_K];
	int flag_occup2[N_K];
	double fac_g1[N_K];
	complex<double> occup[N_VX2];
} *VX_temp;

static int flag_opt_seg;
static int flag_opt_spin;
static int flag_opt_boson;
static int flag_opt_shift;
static int flag_opt_segsp;
static double R_SEG;
static double R_SPIN;
static double R_BOSON;
static double R_SHIFT;
static double R_SEGSP;
static int N_ADD_MIN = 2;
static int N_SHIFT_MIN = 2;

static FILE *fp_log;
static int my_rank=0, process_num=1;
static clock_t time_all_start;

struct wanglandau{
	double log_f;
	double log_gk[N_K];
// 	int kc; // cutoff
} WL;
static const double WL_LOG_F = 1e-10;
static const int WL_KC_MIN = 1;
static const double WL_KC_R = 1.1;
static const int WL_MSR = 2000;
// static const double WL_PK_VAR = 1e-2;


static void sampling(int i_measure, num_mc mc);
static int opt_n_mc(struct num_mc &mc);

static void init_mc_config();
static void init_measure();
static void init_measure_bin();

static void func_measure0();
static void func_measure1();
static void func_measure2();
static void func_measure3();
static void (*func_measure[4])() = {func_measure0, func_measure1, func_measure2, func_measure3};

static void func_averagebin0(int);
static void func_averagebin1(int);
static void func_averagebin2(int);
static void func_averagebin3(int);
static void (*func_averagebin[4])(int) = {func_averagebin0, func_averagebin1, func_averagebin2, func_averagebin3};

static void func_average0(int);
static void func_average1(int);
static void func_average2(int);
static void func_average3(int);
static void (*func_average[4])(int) = {func_average0, func_average1, func_average2, func_average3};

static double mpi_reduce_bin(int);
static double mpi_reduce_accept();


static void add_seg(int sigma, int anti);
static void rem_seg(int sigma, int anti);
// static void (* func_add_rem_seg[2])(int, int)={
// 	add_seg, rem_seg,
// };
static void func_seg(){
	if( rand_int(2) )  add_seg(rand_int(2), rand_int(2));
	else               rem_seg(rand_int(2), rand_int(2));
};

static void add_spin(int anti);
static void rem_spin(int anti);
// static void (* func_spin[2])(int)={
// 	add_spin, rem_spin,
// };
static void func_spin(){
	if( rand_int(2) )  add_spin(rand_int(2));
	else               rem_spin(rand_int(2));
};

static void spin2seg(int dbl_m, int dbl_p);
static void seg2spin(int dbl_m, int dbl_p);
// static void (* func_segspin[2])()={
// 	spin2seg, seg2spin,
// };
static void func_segspin(){
	#if UINF
		if( rand_int(2) )  spin2seg(0, 0);
		else               seg2spin(0, 0);
	#else
		if( rand_int(2) )  spin2seg(rand_int(2), rand_int(2));
		else               seg2spin(rand_int(2), rand_int(2));
	#endif
};

static void perm_boson();
static void shift_tau();
static void change_state0();
static void func_shift(){
	if(X.k)  shift_tau();
	else     change_state0();
};

static void (* func_update[5])() = {
	func_seg, func_spin, perm_boson, func_shift, func_segspin,
};


static void print_config();



//============================================================================
//
// FUNCTIONS TO BE CALLED BY MAIN FUNCTION
//


void sbqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, struct two_particle **p_TP, struct vertex **p_VX)
// void sbqmc_init()
{
	time_all_start = clock();
	
	// A seed of random number (determined from time if seed=0)
	unsigned long seed = RAND_SEED;
	
	#if QMC_MPI
	// 	MPI_Init(&argc_in, &argv_in);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &process_num);
		
		mt_init_by_time_mpi(&seed);
	#else
		mt_init_by_time(&seed);
	#endif // QMC_MPI
	
	if(my_rank==0){
		printf("\nSBQMC_INIT  %s\n", tag);
		printf(" (N_TAU, N_TP, N_TP2) = (%d, %d, %d)\n", N_TAU, N_TP, N_TP2);
		printf(" (N_VX1, N_VX2) = (%d, %d)\n", N_VX1, N_VX2);
		printf(" N_WARMUP = %d\n", N_WARMUP);
		printf(" SELF_VERTEX_ISO = %d\n", SELF_VERTEX_ISO);
		printf(" UINF = %d\n", UINF);
		printf(" seed = %ld\n", seed);
		printf(" IMPROV_ESTIM = %d\n", IMPROV_ESTIM);
		if( VX_NFFT==0 ){
			printf(" VX_NFFT = %d\n", VX_NFFT);
		}
		else{
			printf(" VX_NFFT = %d  M_OS = %d\n", VX_NFFT, M_OS);
		}
	}
	
	#if QMC_MPI
		char host_name[256];
		int host_name_length;
		MPI_Get_processor_name(host_name, &host_name_length);
		
		if(my_rank){
			MPI_Send(host_name, host_name_length+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		}
		else{
			printf("\n MPI parallel with %d nodes\n", process_num);
			printf("  0: %s\n", host_name);
			
			for(int i=1; i<process_num; i++){
				MPI_Status recv_status;
				char other_host_name[256];
				MPI_Recv(other_host_name, 100, MPI_CHAR, i, 0, MPI_COMM_WORLD, &recv_status);
				
				printf(" %2d: %s\n", i, other_host_name);
			}
		}
	#endif
	
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "w");
		fprintf(fp_log, "SB_QMC INIT\n");
		fprintf(fp_log, " mt_init: seed=%ld\n", seed);
		if(QMC_MPI){
			fprintf(fp_log, " MPI parallel with %d nodes\n", process_num);
		}
		fclose(fp_log);
	}
	
	//
	// allocate memory
	//
	for(int s=0; s<2; s++)  G0_alloc(Delta[s], N_TAU);
	D0_alloc(D0_pm, N_TAU);
	for(int a=0; a<2; a++){
		D0_alloc(K0[a], N_TAU);
		D0_alloc(K1[a], N_TAU);
	}
// 	for(int a=0; a<N_A; a++)  S[a] = new cond_op;
// 	
// 	Delta_omega = new complex<double>[N_A][N_TAU/2];
// 	
	B = new phys_quant_bin;
	B_TOT = new phys_quant_bin;
	
	SP = new single_particle [2];
	TP = new two_particle;
	VX = new vertex;
	VX_temp = new vertex_temp[2];
	
	*p_PQ = &PQ;
	*p_SP = SP;
	*p_TP = TP;
	*p_VX = VX;
}
void sbqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, struct two_particle **p_TP)
{
	struct vertex *p_VX;  // dummy
	sbqmc_init(p_PQ, p_SP, p_TP, &p_VX);
}

void sbqmc_set_nmc(struct num_mc n_mc_in)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nSET_NUM_MC\n");
		fclose(fp_log);
	}
	
	n_mc = n_mc_in;
	
	flag_opt_seg = 0;
	flag_opt_spin = 0;
	flag_opt_boson = 0;
	flag_opt_shift = 0;
	flag_opt_segsp = 0;
	
	if(n_mc.seg<0){
		flag_opt_seg = 1;
		R_SEG = - (double)n_mc.seg / 10.;
	}
	if(n_mc.spin<0){
		flag_opt_spin = 1;
		R_SPIN = - (double)n_mc.spin / 10.;
	}
	if(n_mc.boson<0){
		flag_opt_boson = 1;
		R_BOSON = - (double)n_mc.boson / 10.;
	}
	if(n_mc.shift<0){
		flag_opt_shift = 1;
		R_SHIFT = - (double)n_mc.shift / 10.;
	}
	if(n_mc.segsp<0){
		flag_opt_segsp = 1;
		R_SEGSP = - (double)n_mc.segsp / 10.;
	}
}

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
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double *fac_int, double *fac_chi_in, double D_jump)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nSET_INPUT\n");
		fclose(fp_log);
	}
	
	prm = prm_in;
	
// 	if( V_sqr[0] == 0 && V_sqr[1] == 0 ){
	if( Delta_omega_in == NULL ){
		flag_V = 0;
		for(int s=0; s<2; s++){
			for(int i=0; i<N_TAU/2; i++)  Delta_omega[s][i] = 0;
		}
	}
	else{
		flag_V = 1;
		for(int s=0; s<2; s++){
			for(int i=0; i<N_TAU/2; i++)  Delta_omega[s][i] = Delta_omega_in[s][i];
		}
	}
// 	if( LIMIT_PARAM != 2 ){  // V != 0
	if( flag_V ){
		for(int s=0; s<2; s++){
// 			for(int i=0; i<N_TAU/2; i++)  Delta_omega[s][i] = Delta_omega_in[s][i];
			G0_init_fft(Delta[s], Delta_omega[s], prm.beta, V_sqr[s]);
		}
	}
	
	for(int a=0; a<3; a++){
		if( D0_omega_in[a] == NULL ){
			flag_g[a] = 0;
			for(int i=0; i<=N_TAU/2; i++)  D0_omega[a][i] = 0;
		}
		else{
			double fac[3] = {1, 4, 4};  // Pauli spin to 1/2-spin
			flag_g[a] = 1;
// 			for(int i=0; i<=N_TAU/2; i++)  D0_omega[a][i] = D0_omega_in[a][i];
			for(int i=0; i<=N_TAU/2; i++)  D0_omega[a][i] = D0_omega_in[a][i] * fac[a] * fac_int[a];
		}
		fac_chi[a] = fac_chi_in[a];
	}
	
	// ch, sp(zz)
	for(int a=0; a<2; a++){
		if( flag_g[a] ){
			K0_init_sum(K0[a], D0_omega[a], prm.beta);
			K1_init_sum(K1[a], D0_omega[a], prm.beta);
		}
	}
	
	for(int s=0; s<2; s++){
		prm_ef_unrenorm[s] = prm.ef[s];
		prm.ef[s] -= real(D0_omega[0][0]) / 2.;  // ch
		prm.ef[s] += real(D0_omega[1][0]) / 8.;  // sp
	}
	prm_U_unrenorm = prm.U;
	prm.U += real(D0_omega[0][0]);       // ch
	prm.U -= real(D0_omega[1][0]) / 4.;  // sp
	
	// sp(pm)
// 	if( LIMIT_PARAM != 1 ){  // g_xy != 0
	if( flag_g[2] ){
// 		D0_init_fft(D0_pm, D0_omega[2], prm.beta, D_jump);
		D0_init_fft(D0_pm, D0_omega[2], prm.beta, D_jump * 4. * fac_int[2]);
	}
	
	if(my_rank==0 && DISPLAY){
		printf("\nSBQMC_SET_INPUT\n");
		printf(" fac_int = %.3lf, %.3lf, %.3lf\n", fac_int[0], fac_int[1], fac_int[2]);
		printf(" fac_chi = %.3lf, %.3lf, %.3lf\n", fac_chi[0], fac_chi[1], fac_chi[2]);
		printf(" ef = %.4e, %.4e\n", prm.ef[0], prm.ef[1]);
		printf(" U  = %.4e\n", prm.U);
		printf(" flag(V, g_ch, g_sp, g_pm) = (%d, %d, %d, %d)\n",
		 flag_V, flag_g[0], flag_g[1], flag_g[2]);
	}
}
// D_jump = 0
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double *fac_int, double *fac_chi_in)
{
	sbqmc_set_input(prm_in, Delta_omega_in, V_sqr, D0_omega_in, fac_int, fac_chi_in, 0);
}
void sbqmc_set_input(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in)
{
	double fac_int[3] = {1,1,1};
	double fac_chi[3] = {1,1,1};
	sbqmc_set_input(prm_in, Delta_omega_in, V_sqr, D0_omega_in, fac_int, fac_chi, 0);
}

/*
//  Delta = V_sqr * G0
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> G0_omega[2][N_TAU/2], double V_sqr[2], complex<double> **D0_omega_in, double D_jump)
{
	complex<double> (*Delta_omega_in)[N_TAU/2] = new complex<double>[2][N_TAU/2];
	
	for(int a=0; a<2; a++){
		for(int i=0; i<N_TAU/2; i++)  Delta_omega_in[a][i] = G0_omega[a][i] * V_sqr[a];
	}
	
	sbqmc_set_input1(prm_in, Delta_omega_in, V_sqr, D0_omega_in, D_jump);
	
	delete [] Delta_omega_in;
}
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> G0_omega[2][N_TAU/2], double V_sqr[2], complex<double> **D0_omega_in)
{
	sbqmc_set_input2(prm_in, G0_omega, V_sqr, D0_omega_in, 0);
}
*/

/*
// interactions are defined with the Pauli operator:  g * sigma . phi
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in, double D_jump)
{
	if(my_rank==0 && DISPLAY){
		printf("\nINPUT2\n");
		printf(" D0 --> 4*D0  (sigma.phi to S.phi)\n");
	}
	
	complex<double> (*D0_omega_z) = new complex<double>[N_TAU/2+1];
	complex<double> (*D0_omega_xy) = new complex<double>[N_TAU/2+1];
	
	complex<double> *D0_omega_in2[3] = {D0_omega_in[0], D0_omega_z, D0_omega_xy};
	
	// factor 4: Pauli operator to 1/2-spin
	for(int a=1; a<3; a++){
		if( D0_omega_in[a] == NULL ){
			D0_omega_in2[a] = NULL;
		}
		else{
			for(int i=0; i<=N_TAU/2; i++){
				D0_omega_in2[a][i] = D0_omega_in[a][i] * 4.;
			}
		}
	}
	D_jump *= 4.;
	
	sbqmc_set_input1(prm_in, Delta_omega_in, V_sqr, D0_omega_in2, D_jump);
	
	delete [] D0_omega_z;
	delete [] D0_omega_xy;
}
void sbqmc_set_input2(struct sb_qmc_params prm_in, complex<double> **Delta_omega_in, double V_sqr[2], complex<double> **D0_omega_in)
{
	sbqmc_set_input2(prm_in, Delta_omega_in, V_sqr, D0_omega_in, 0);
}
*/

// 0: successful exit
// 1: ave_sign is lower than MIN_AVE_SIGN
// 2: average expansion order is larger than MAX_AVE_PERT
int sbqmc_eval(int flag_tp)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nEVAL_MC\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nSBQMC_EVAL\n");
	}
	time_t time_start = clock();
	
	init_mc_config();
	init_measure();
	
	//
	// warming up
	//
	struct num_mc n_warmup;
	n_warmup.msr = N_WARMUP * process_num;
	n_warmup.bin = 1;
	n_warmup.seg = 1;
	n_warmup.spin = 1;
	n_warmup.boson = 1;
	n_warmup.shift = 1;
	n_warmup.segsp = 1;
	
	if( n_mc.shift == 0 )  n_warmup.shift = 0;
	if( n_mc.segsp == 0 )  n_warmup.segsp = 0;
	
	if( ! flag_g[2] )  n_warmup.spin = n_warmup.boson = n_warmup.segsp = 0;  // g_xy=0
	if( ! flag_V    )  n_warmup.seg = n_warmup.segsp = 0;  // V=0
	
// 	if(R_SHIFT)  sampling(0, 1, N_WARMUP * process_num, 1, 1, 1);
// 	else         sampling(0, 1, N_WARMUP * process_num, 1, 1, 0);
	sampling(0, n_warmup);
	
// 	if( PQ.ave_Z_V[0] > MAX_AVE_PERT_V || PQ.ave_Z_V[1] > MAX_AVE_PERT_V )  return(2);
// 	if( PQ.ave_Z_g > MAX_AVE_PERT_G )  return(2);
	
	//
	// optimizing parameters
	//
// 	if( flag_opt_add || flag_opt_exc || flag_opt_shift ){
	{
		int flag_return = 0;
		if(my_rank==0)  flag_return = opt_n_mc(n_mc);
		#if QMC_MPI
			MPI_Bcast(&flag_return, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&n_mc, 7, MPI_INT, 0, MPI_COMM_WORLD);
		#endif
		if( flag_return ){
			if(my_rank==0 && DISPLAY)  printf(" *** return %d\n", flag_return);
			return( flag_return );
		}
	}
	if( ! flag_g[2] )  n_mc.spin = n_mc.boson = n_mc.segsp = 0;  // g=0
	if( ! flag_V    )  n_mc.seg = n_mc.segsp = 0;  // V=0
	
	//
	// measuring physical quantities
	//
	init_measure();
	
	sampling(flag_tp+1, n_mc);
	
	if(my_rank==0 && DISPLAY){
		time_t time_end = clock();
		char str[100];
		sprint_time(str, time_end - time_start);
		printf("\n time:%s", str);
	}
	
	return(0);
}

// evaluate free energy by Wang-Landau sampling
int sbqmc_eval_thermo()
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nEVAL_THERMO\n");
		fclose(fp_log);
	}
	
// 	if( LIMIT_PARAM != 2){
	if( flag_V ){
		printf("*** NOT IMPLEMENTED\n");
		exit(0);
	}
	
	if(my_rank==0 && DISPLAY){
		printf("\n WL_LOG_F = %.2e\n", WL_LOG_F);
		printf(" WL_KC_MIN = %d\n", WL_KC_MIN);
		printf(" WL_KC_R = %.3lf\n", WL_KC_R);
		printf(" WL_MSR = %d\n", WL_MSR);
	}
	
	init_mc_config();
	init_measure();
	
	//
	// warming up
	//
	struct num_mc n_warmup;
	n_warmup.msr = N_WARMUP * process_num;
	n_warmup.bin = 1;
	n_warmup.seg = 1;
	n_warmup.spin = 1;
	n_warmup.boson = 1;
	n_warmup.shift = 1;
	n_warmup.segsp = 1;
	
	if( n_mc.shift == 0 )  n_warmup.shift = 0;
	if( n_mc.segsp == 0 )  n_warmup.segsp = 0;
	
	if( ! flag_g[2] )  n_warmup.spin = n_warmup.boson = n_warmup.segsp = 0;  // g=0
	if( ! flag_V    )  n_warmup.seg = n_warmup.segsp = 0;  // V=0
	
	sampling(0, n_warmup);
	
	if( PQ.ave_Z_V[0] > MAX_AVE_PERT_V || PQ.ave_Z_V[1] > MAX_AVE_PERT_V )  return(2);
	if( PQ.ave_Z_g > MAX_AVE_PERT_G )  return(2);
	
	//
	// optimizing parameters
	//
// 	if( flag_opt_add || flag_opt_exc || flag_opt_shift ){
	{
		int flag_return = 0;
		if(my_rank==0)  flag_return = opt_n_mc(n_mc);
		#if QMC_MPI
			MPI_Bcast(&flag_return, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&n_mc, 7, MPI_INT, 0, MPI_COMM_WORLD);
		#endif
		if( flag_return ){
			if(my_rank==0 && DISPLAY)  printf(" *** return %d\n", flag_return);
			return(2);
		}
	}
	if( ! flag_g[2] )  n_mc.spin = n_mc.boson = n_mc.segsp = 0;  // g=0
	if( ! flag_V    )  n_mc.seg = n_mc.segsp = 0;  // V=0
	
	//
	// tuning WL factor WL.log_gk[]
	//
	
	// init
	init_mc_config();
// 	for(int i=0; i<N_K; i++)  WL.log_gk[i] = 0;
	
	int WL_kc = (int)(PQ.ave_Z_g * WL_KC_R);
	if( WL_kc < WL_KC_MIN )  WL_kc = WL_KC_MIN;
	if(my_rank==0){
		printf("\n WL_kc = %d\n", WL_kc);
	}
	#if QMC_MPI
		MPI_Bcast(&WL_kc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif
	for(int i=WL_kc+1; i<N_K; i++)  WL.log_gk[i] = DBL_MAX;
	
	n_warmup.msr = (n_mc.spin + n_mc.boson) * WL_MSR;
	
	WL.log_f = 1.0;
	do{
		if(my_rank==0){
			printf("\n--- log_f = %.4e\n", WL.log_f);
		}
		
		init_measure();
		sampling(0, n_warmup);
		
		for(int i=0; i<=WL_kc; i++)  WL.log_gk[WL_kc-i] -= WL.log_gk[0];
		
		#if QMC_MPI
			double temp[WL_kc+1];
			MPI_Allreduce(WL.log_gk, temp, WL_kc+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(int i=0; i<=WL_kc; i++)  WL.log_gk[i] = temp[i] / (double)process_num;
		#endif
		
// 		printf("\n log_gk[]\n");
// 		for(int i=0; i<=WL_kc; i++){
// 			printf("  %d %lf\n", i, WL.log_gk[i]);
// 		}
		
		// check if the histogram is flat
		if(my_rank==0){
// 			double max_pk = PQ.Z_g[0];
// 			double min_pk = PQ.Z_g[0];
// 			for(int i=1; i<=WL_kc; i++){
// 				if( max_pk < PQ.Z_g[i] )  max_pk = PQ.Z_g[i];
// 				if( min_pk > PQ.Z_g[i] )  min_pk = PQ.Z_g[i];
// 			}
// 			printf("\n measure of flatness: %lf\n", min_pk / max_pk);
			
			double pk_var = 0, pk_mean = 0;
			for(int i=0; i<=WL_kc; i++){
				pk_mean += PQ.Z_g[i] / double(WL_kc+1);
				pk_var += pow(PQ.Z_g[i], 2) / double(WL_kc+1);
			}
			pk_var = sqrt( pk_var - pow(pk_mean, 2) );
			printf("\n measure of flatness: %.1e\n", pk_var);
			
// 			if( pk_var > WL_PK_VAR ){
// 				n_warmup.msr *= 2;
// 				printf("  msr *= 2\n");
// // 				WL.log_f *= 2.;
// 			}
		}
		#if QMC_MPI
			MPI_Bcast(&n_warmup.msr, 1, MPI_INT, 0, MPI_COMM_WORLD);
		#endif
		
		WL.log_f *= 0.5;
	} while( WL.log_f > WL_LOG_F );
	
	WL.log_f = 0;
	for(int i=WL_kc+1; i<N_K; i++)  WL.log_gk[i] = WL.log_gk[WL_kc];
	
	//
	// measuring physical quantities
	//
	init_measure();
	sampling(0, n_mc);
	
// 	printf("\n free_energy = %.5e +- %.4e\n", PQ.free_energy, PQ.free_energy_err);
	
	return(0);
}

void sbqmc_final()
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nFINAL\n");
		fclose(fp_log);
	}
	
	if(my_rank==0 && DISPLAY){
		printf("\nSBQMC_FINAL\n");
		char str[100];
		sprint_time_mpi(str, time_all_start, clock(), 0.0);
		printf("%s", str);
		sbqmc_fprint_log(str);
	}
	
	
	//
	// free memory
	//
	for(int s=0; s<2; s++)  G0_free(Delta[s]);
	D0_free(D0_pm);
	for(int a=0; a<2; a++){
		D0_free(K0[a]);
		D0_free(K1[a]);
	}
// 	for(int a=0; a<N_A; a++)  delete S[a];
// 	delete [] Delta_omega;
	
	delete B;
	delete B_TOT;
	
	delete [] SP;
	delete TP;
	delete VX;
	delete [] VX_temp;
}

void sbqmc_fprint_log(char *str)
{
	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "%s", str);
	fclose(fp_log);
}


//============================================================================
//
// MC SAMPLING
//

static void eval_acceptance_sub(double tot, double n_accept, double n_reject, double &d_accept, double &d_reject, char *str)
{
// 	double tot_seg = n_sample * (double)n_seg;
	d_accept = n_accept / tot;
	d_reject = n_reject / tot;
	
	if(DISPLAY){
		printf("%s", str);
		printf(" accept %.6lf | reject %.6lf\n", d_accept, d_reject);
	}
	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "%s", str);
	fprintf(fp_log, " accept %.6lf | reject %.6lf\n", d_accept, d_reject);
	fclose(fp_log);
}

static void eval_acceptance(double n_sample, int n_seg, int n_spin, int n_boson, int n_shift, int n_segsp)
{
	//
	// acceptance
	//
	char str[32];
	
	sprintf(str, "\n add/rem segment :");  // 18 letters
	eval_acceptance_sub(n_sample * (double)n_seg, (double)ACCPT_sub.n_accept_seg, (double)ACCPT_sub.n_reject_seg,
	 ACCPT.d_accept_seg, ACCPT.d_reject_seg, str);
	
	sprintf(str, " add/rem spin    :");
	eval_acceptance_sub(n_sample * (double)n_spin, (double)ACCPT_sub.n_accept_spin, (double)ACCPT_sub.n_reject_spin,
	 ACCPT.d_accept_spin, ACCPT.d_reject_spin, str);
	
	sprintf(str, " permanent boson :");
	eval_acceptance_sub(n_sample * (double)n_boson, (double)ACCPT_sub.n_accept_boson, (double)ACCPT_sub.n_reject_boson,
	 ACCPT.d_accept_boson, ACCPT.d_reject_boson, str);
	
	if(n_shift){
		sprintf(str, " shift           :");
		eval_acceptance_sub(n_sample * (double)n_shift, (double)ACCPT_sub.n_accept_shift, (double)ACCPT_sub.n_reject_shift,
		 ACCPT.d_accept_shift, ACCPT.d_reject_shift, str);
	}
	else{
		ACCPT.d_accept_shift = 0;
		ACCPT.d_reject_shift = 0;
	}
	
	if(n_segsp){
		sprintf(str, " change seg&spin :");
		eval_acceptance_sub(n_sample * (double)n_segsp, (double)ACCPT_sub.n_accept_segsp, (double)ACCPT_sub.n_reject_segsp,
		 ACCPT.d_accept_segsp, ACCPT.d_reject_segsp, str);
	}
	else{
		ACCPT.d_accept_segsp = 0;
		ACCPT.d_reject_segsp = 0;
	}
}

static void sampling(int i_measure, num_mc mc)
{
// 	unsigned long local_n_sample = (n_sample + process_num - 1)/(process_num);  // rounding upward
	unsigned long local_n_sample = (mc.msr + process_num - 1)/(process_num);  // rounding upward
	
	unsigned long global_n_sample = local_n_sample * process_num;
	
	if(my_rank==0){
		if(DISPLAY)  printf("\nSAMPLING\n");
		
		char str[128];
		sprintf(str, "\n n_bin   : %d\n", mc.bin);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_sample: %ld  (local:%ld)\n", global_n_sample, local_n_sample);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_seg   : %d\n", mc.seg);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_spin  : %d\n", mc.spin);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_boson : %d\n", mc.boson);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_shift : %d\n", mc.shift);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
		
		sprintf(str, " n_segsp : %d\n", mc.segsp);
		if(DISPLAY)  printf("%s", str);
		sbqmc_fprint_log(str);
	}
	
	//
	// init acceptance
	//
	ACCPT_sub.n_accept_seg = ACCPT_sub.n_reject_seg = 0;
	ACCPT_sub.n_accept_spin = ACCPT_sub.n_reject_spin = 0;
	ACCPT_sub.n_accept_shift = ACCPT_sub.n_reject_shift = 0;
	ACCPT_sub.n_accept_boson = ACCPT_sub.n_reject_boson = 0;
	ACCPT_sub.n_accept_segsp = ACCPT_sub.n_reject_segsp = 0;
	
	
/*	
// 	print_config();
	
	for(int i=0; i<5; i++){
// 		add_seg(0);
		add_spin();
	}
	for(int i=0; i<5; i++){
// 		add_seg(0);
// 		rem_seg(0);
// 		add_spin();
		rem_spin();
	}
// 	for(int i=0; i<100; i++){
// 		func_add_rem_seg[rand_int(2)](0);
// 	}
// 	
	print_config();
	exit(0);
*/	
	
	
	//
	// sampling
	//
	//static void (* func_update[5])() = {
	//	func_seg, func_spin, perm_boson, func_shift, func_segspin,
	//};
	int n_update = mc.seg + mc.spin + mc.boson + mc.shift + mc.segsp;
	int update_type[n_update];
	{
		int n = 0;
		for(int l=0; l<mc.seg; l++)    update_type[n++] = 0;
		for(int l=0; l<mc.spin; l++)   update_type[n++] = 1;
		for(int l=0; l<mc.boson; l++)  update_type[n++] = 2;
		for(int l=0; l<mc.shift; l++)  update_type[n++] = 3;
		for(int l=0; l<mc.segsp; l++)  update_type[n++] = 4;
	}
	
	if(my_rank==0 && DISPLAY){
		printf("\n|---------|---------|---------|---------|---------|\n|");
	}
	unsigned long n_meter = (unsigned long)(mc.bin * local_n_sample / 50);
	unsigned long i_meter = 0;
	
	clock_t time_start = clock();
	double time_trans_tot=0;
	
	for(int i=0; i<mc.bin; i++){
		clock_t time_bin_start = clock();
		double time_trans_bin=0;
		
		init_measure_bin();
		
		for(unsigned long j=0; j<local_n_sample; j++){
			
			int perm[n_update];
			random_permutation(perm, n_update);
			for(int l=0; l<n_update; l++)  func_update[ update_type[ perm[l] ] ]();
			
// 			for(int i_seg=0; i_seg<mc.seg; i_seg++){
// 				func_add_rem_seg[rand_int(2)](rand_int(2), rand_int(2));
// // 				func_add_rem_seg[rand_int(2)](0, rand_int(2));  // spinless
// 			}
// 			for(int i_spin=0; i_spin<mc.spin; i_spin++){
// 				func_spin[rand_int(2)](rand_int(2));
// 				for(int i_boson=0; i_boson<mc.boson; i_boson++)  perm_boson();
// 			}
// 			for(int i_segsp=0; i_segsp<mc.segsp; i_segsp++){
// 				func_segspin[rand_int(2)]();
// 			}
// 			if(X.k){
// 				for(int l=0; l<mc.shift; l++){
// 					shift_tau();
// 				}
// 			} else {
// 				change_state0();
// 			}
			
			func_measure[i_measure]();
			
			i_meter++;
			if(i_meter == n_meter){
				i_meter = 0;
				if(my_rank==0 && DISPLAY){
					printf("*");
					fflush(stdout);
				}
			}
		}
		
		
		#if QMC_MPI
		time_trans_bin = mpi_reduce_bin(i_measure);
		time_trans_tot += time_trans_bin;
		#endif // QMC_MPI
		
		
		if(my_rank==0){
			func_averagebin[i_measure](global_n_sample);
			
			clock_t time_bin_end = clock();
	// 		print_time(time_start, time_end);
// 			print_time(time_bin_start, time_bin_end, LOG_FILE);
// 			print_time_mpi(time_bin_start, time_bin_end, time_trans_bin, LOG_FILE);
			
// 			char str[100];
// 			sprint_time_mpi(str, time_bin_start, time_bin_end, time_trans_bin);
// 			sbqmc_fprint_log(str);
		}
	}
	if(my_rank==0 && DISPLAY){
		printf("\n");
		fflush(stdout);
	}
	
	#if QMC_MPI
	time_trans_tot += mpi_reduce_accept();
	#endif // QMC_MPI
	
	if(my_rank==0){
		eval_acceptance(global_n_sample * mc.bin, mc.seg, mc.spin, mc.boson, mc.shift, mc.segsp);
		func_average[i_measure](mc.bin);
	}
}

// return 2 if the average expansion order exceeds the limit
static int opt_n_mc(struct num_mc &mc)
{
	if( PQ.ave_Z_V[0] > MAX_AVE_PERT_V )  return(2);
	if( PQ.ave_Z_V[1] > MAX_AVE_PERT_V )  return(2);
	if( PQ.ave_Z_g > MAX_AVE_PERT_G )  return(2);
	
	//
	// determine mc.??? if( flag_opt_??? )
	//
	if( flag_opt_seg ){
		double ave_Z_V = PQ.ave_Z_V[0] + PQ.ave_Z_V[1];
		mc.seg = (int)(ave_Z_V / ACCPT.d_accept_seg * R_SEG);
// 		printf(" mc.seg=%d\n", mc.seg);
		if(mc.seg < N_ADD_MIN)  mc.seg = N_ADD_MIN;
	}
	
	if( flag_opt_spin ){
		mc.spin = (int)(PQ.ave_Z_g / ACCPT.d_accept_spin * R_SPIN);
// 		printf(" mc.spin=%d\n", mc.spin);
		if(mc.spin < N_ADD_MIN)  mc.spin = N_ADD_MIN;
	}
	
	if( flag_opt_boson ){
// 		mc.boson = (int)(PQ.ave_Z_g / ACCPT.d_accept_boson * R_BOSON);
// 		mc.boson = (int)(PQ.ave_Z_g * R_BOSON);
		mc.boson = (int)(PQ.ave_Z_g * R_BOSON * (double)mc.spin);
// 		printf(" mc.boson=%d\n", mc.boson);
		if(mc.boson < N_ADD_MIN)  mc.boson = N_ADD_MIN;
	}
	
	if( flag_opt_shift ){
		double ave_tot_k = PQ.ave_Z_V[0] + PQ.ave_Z_V[1] + PQ.ave_Z_g;
		mc.shift = (int)(ave_tot_k / ACCPT.d_accept_shift * R_SHIFT);
// 		printf(" mc.shift=%d\n", mc.shift);
		if(mc.shift < N_SHIFT_MIN)  mc.shift = N_SHIFT_MIN;
	}
	
	if( flag_opt_segsp ){
// 		double min_ave_Z = MIN( 0.5 * (PQ.ave_Z_V[0] + PQ.ave_Z_V[1]), PQ.ave_Z_g );
		// geometrical average
		double ave_Z = sqrt( 0.5 * (PQ.ave_Z_V[0] + PQ.ave_Z_V[1]) *  PQ.ave_Z_g );
		mc.segsp = (int)(ave_Z * R_SEGSP);
// 		printf(" mc.segsp=%d\n", mc.segsp);
		if(mc.segsp < N_ADD_MIN)  mc.segsp = N_ADD_MIN;
	}
	return 0;
}

//============================================================================
//
// INITIALIZE FOR MC SAMPLING
//

static void init_mc_config()
{
	w_sign = +1;

	X.k = 0;
// 	X.op[0] = F_UP_CR;  // op2state[F_UP_CR] == EMP
// 	X.op[0] = F_UP_AN; // op2state[F_UP_AN] == UP
// 	X.op[0] = S_PLUS;  // op2state[S_PLUS] == DW
// 	X.tau[0] = 0;
// 	X.tau[1] = prm.beta;
// 	X.tau[0] = prm.beta;
	X.state[0] = UP;
	
	S[0].k = S[1].k = 0;
	Y.k = 0;
	
	
	WL.log_f = 0;
	for(int i=0; i<N_K; i++)  WL.log_gk[i] = 0;
	
	
	double delta_tau2 = prm.beta / (double)(2*N_TP2);
	
	double TP2_tau[N_TP2+1];
	for(int i=0; i<=N_TP2; i++)  TP2_tau[i] = delta_tau2 * (double)i;
	
	TP->tau[0] = TP2_tau[0];
	int R_TP = N_TP2/N_TP;
	int i1=0, i2=0;
	for(int j=0; j<R_TP; j++){
		for(int i=0; i<N_TP/2/R_TP; i++){
			i1 ++;
			i2 += j+1;
			TP->tau[i1] = TP2_tau[i2];
		}
	}
	for(int j=0; j<R_TP; j++){
		for(int i=0; i<N_TP/2/R_TP; i++){
			i1 ++;
			i2 += j+R_TP;
			TP->tau[i1] = TP2_tau[i2];
		}
	}
}

static void init_measure()
{
// 	struct phys_quant{
// 		double ave_sign, ave_sign_err;
// 		double Z_V[2][N_K], Z_V_err[2][N_K];  // V
// 		double Z_g[N_K], Z_g_err[N_K];  // g
// 		double ave_Z_V[2], ave_Z_V_err[2];
// 		double ave_Z_g, ave_Z_g_err;
// 		double occup_tot, occup_tot_err;
// 		double occup_dif, occup_dif_err;
// 		double occup_ud[2], occup_ud_err[2];
// 		double occup_dbl, occup_dbl_err;
// 		double stat_suscep_sp, stat_suscep_sp_err;
// 		double stat_suscep_ch, stat_suscep_ch_err;
// 		double stat_suscep_pm, stat_suscep_pm_err;
// 	};
	
	PQ.ave_sign = 0;
	PQ.ave_sign_err = 0;
	
	for(int s=0; s<2; s++){
		for(int i=0; i<N_K; i++)  PQ.Z_V[s][i] = PQ.Z_V_err[s][i] = 0;
	}
	for(int i=0; i<N_K; i++)  PQ.Z_g[i] = PQ.Z_g_err[i] = 0;
	
	for(int s=0; s<2; s++){
		PQ.ave_Z_V[s] = PQ.ave_Z_V_err[s] = 0;
	}
	PQ.ave_Z_g = PQ.ave_Z_g_err = 0;
	
	PQ.occup_tot = PQ.occup_tot_err = 0;
	PQ.occup_dif = PQ.occup_dif_err = 0;
	for(int s=0; s<2; s++){
		PQ.occup_ud[s] = PQ.occup_ud_err[s] = 0;
	}
	PQ.occup_dbl = PQ.occup_dbl_err = 0;
	
	PQ.stat_suscep_sp = PQ.stat_suscep_sp_err = 0;
	PQ.stat_suscep_ch = PQ.stat_suscep_ch_err = 0;
	PQ.stat_suscep_pm = PQ.stat_suscep_pm_err = 0;
	
	PQ.m_z_sqr = PQ.m_z_sqr_err = 0;
	PQ.m_z_pw4 = PQ.m_z_pw4_err = 0;
	PQ.m_xy_sqr = PQ.m_xy_sqr_err = 0;
	PQ.m_xy_pw4 = PQ.m_xy_pw4_err = 0;
	PQ.m_rt_sqr = PQ.m_rt_sqr_err = 0;
	PQ.m_rt_pw4 = PQ.m_rt_pw4_err = 0;
	
	PQ.free_energy = PQ.free_energy_err = 0;
	
	//
	// SP
	//
	for(int s=0; s<2; s++){
		for(int i=0; i<=N_TAU; i++){
			SP[s].Gf_tau[i] = SP[s].Gf_tau_err[i] = 0;
			SP[s].G_self[i] = 0;
		}
	}
	
	//
	// TP
	//
	for(int i=0; i<=N_TP; i++){
		TP->chi_sp_tau[i] = TP->chi_sp_tau_err[i] = 0;
		TP->chi_ch_tau[i] = TP->chi_ch_tau_err[i] = 0;
	}
	
	for(int i=0; i<=N_TAU; i++){
		TP->chi_pm_tau[i] = TP->chi_pm_tau_err[i] = 0;
	}
	
	//
	// VX
	//
// 	for(int iw1=0; iw1<2*N_VX1; iw1++){
// 		for(int iw2=0; iw2<2*N_VX1; iw2++){
// 			for(int nu=0; nu<N_VX2; nu++){
// 				VX->chi_ch[iw1][iw2][nu] = 0;
// 				VX->chi_sp[iw1][iw2][nu] = 0;
// 				VX->chi_pm[iw1][iw2][nu] = 0;
// 			}
// 		}
// 	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						VX->Gfour_lo[s1][s2][iw1][iw2][nu] = 0;
						VX->gamma_lo[s1][s2][iw1][iw2][nu] = 0;
					}
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					VX->Gfour_tr[s][iw1][iw2][nu] = 0;
					VX->gamma_tr[s][iw1][iw2][nu] = 0;
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			VX->G[s][iw] = VX->self[s][iw] = 0;
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int nu=0; nu<N_VX2; nu++){
				VX->suscep_lo[s1][s2][nu] = 0;
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					VX->Gthree_lo[s1][s2][iw][nu] = 0;
					VX->lambda_lo[s1][s2][iw][nu] = 0;
				}
			}
		}
	}
}

static void init_measure_bin()
{
// 	static struct phys_quant_bin{
// 		unsigned long n_k[2][N_K];  // V
// 		unsigned long n_l[N_K];  // g
// 		long int ave_sign;
// 		double Gf[2][N_TAU+1];
// 		double occup[4];
// 		double chi[4][4][N_TP+1];
// 		double chi_pm[N_TAU+1];
// 	// 	double chi_tr1[2*N_TP2+1], chi_tr2[2*N_TP2+1];
// 	} *B, *B_TOT;
	
	B->ave_sign = 0;
	
	for(int i=0; i<N_K; i++){
		for(int s=0; s<2; s++)  B->n_k[s][i] = 0;
		B->n_l[i] = 0;
	}
	
	B->occup[0] = B->occup[1] = B->occup[2] = B->occup[3] = 0;
	for(int s=0; s<2; s++){
		for(int i=0; i<=N_TAU; i++){
			B->Gf[s][i] = 0;
			B->self[s][i] = 0;
		}
	}
	
	B->m_z_sqr = B->m_z_pw4 = 0;
	B->m_xy_sqr = B->m_xy_pw4 = 0;
	B->m_rt_sqr = B->m_rt_pw4 = 0;
	
	for(int i=0; i<=N_TAU; i++)  B->chi_pm[i] = 0;
	
	
	for(int a1=0; a1<4; a1++){
		for(int a2=0; a2<4; a2++){
			for(int i=0; i<=N_TP; i++)  B->chi[a1][a2][i] = 0;
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						B->vx_G4_lo[s1][s2][iw1][iw2][nu] = 0;
						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] = 0;
					}
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G4_tr[s][iw1][iw2][nu] = 0;
					B->vx_gamma_tr[s][iw1][iw2][nu] = 0;
				}
			}
		}
	}
	
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			B->vx_Gf[s][iw] = B->vx_self[s][iw] = 0;
		}
	}
	
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int nu=0; nu<N_VX2; nu++){
				B->vx_suscep[s1][s2][nu] = 0;
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G3[s1][s2][iw][nu] = 0;
				}
			}
		}
	}
}

//============================================================================
//
// MEASUREMENT
//

// return 1 if spin sigma is occupied
// return 0 if               unoccupied
static int occupied(int state, int sigma)
{
	if( state & (0x01 << sigma) )  return 1;
	else  return 0;
}

// return +1 if s == 0
//        -1    s == 1
static int sign_spin1(int s)
{
	return (-2*s+1);
}
// return +1 if s1 == s2
//        -1    s1 != s2
static int sign_spin2(int s1, int s2)
{
	return ( -2 * abs(s1-s2) + 1 );
}

static inline void measure_stat()
{
	for(int s=0; s<2; s++){
		B->n_k[s][S[s].k] ++;
	}
	B->n_l[Y.k] ++;
	
	B->ave_sign += w_sign;
	
	WL.log_gk[Y.k] += WL.log_f;
}

static double mz_sqr;  // calculated in measure_sp() and used in measure_pm()

static inline void measure_sp()
{
	//
	// occupation number
	//
	double temp_occup[4] = {0};
	if(X.k){
		for(int i=0; i<X.k; i++){
// 			B->occup[op2state[X.op[i+1]]] += (X.tau[i+1] - X.tau[i]) * (double)w_sign;
// 			temp_occup[op2state[X.op[i+1]]] += (X.tau[i+1] - X.tau[i]);
			temp_occup[X.state[i+1]] += (X.tau[i+1] - X.tau[i]);
		}
	}
	else{
// 		B->occup[op2state[X.op[0]]] += prm.beta * (double)w_sign;
// 		temp_occup[op2state[X.op[0]]] += prm.beta;
		temp_occup[X.state[0]] += prm.beta;
	}
	for(int a=0; a<4; a++)  B->occup[a] += temp_occup[a] * (double)w_sign;
	VX_temp[0].occup[0] = temp_occup[UP] + temp_occup[DBL];
	VX_temp[1].occup[0] = temp_occup[DW] + temp_occup[DBL];
	
	double mz = temp_occup[UP] - temp_occup[DW];
// 	B->mom = mom * (double)w_sign;
	mz_sqr = mz * mz / (prm.beta * prm.beta);
	B->m_z_sqr += mz_sqr * (double)w_sign;
	B->m_z_pw4 += mz_sqr * mz_sqr * (double)w_sign;
	
	//
	// Green function
	//
	double delta_tau = prm.beta / (double)N_TAU;
	
	for(int s=0; s<2; s++){
		// for self-energy (from static U)
// 		int flag_occup1[S[s].k], flag_occup2[S[s].k];
		for(int i=0; i<S[s].k; i++){
			int i_tau1 = tau_position(X.tau, X.k, S[s].tau1[i]);
			int i_tau2 = tau_position(X.tau, X.k, S[s].tau2[i]);
			VX_temp[s].flag_occup1[i] = occupied(X.state[i_tau1], (s+1)%2) ? 1 : 0;
			VX_temp[s].flag_occup2[i] = occupied(X.state[i_tau2], (s+1)%2) ? 1 : 0;
		}
		
		// for self-energy (from g_ch, g_sp)
// 		double fac_g1[S[s].k];
		double fac_init = flag_g[0] != 0 ? -real(D0_omega[0][0]) : 0;
		for(int i=0; i<S[s].k; i++)  VX_temp[s].fac_g1[i] = fac_init;
		
		double spin_fac = SELF_VERTEX_ISO ? 3./4. : 1./4.;
		for(int a=0; a<2; a++){
			if( flag_g[a] ){
				for(int i=0; i<S[s].k; i++){
					for(int ss=0; ss<2; ss++){
						double temp = 0;
						for(int j=0; j<S[ss].k; j++){
							temp -= D0_calc_interp(K1[a], S[s].tau1[i], S[ss].tau1[j]);  // f-
							temp += D0_calc_interp(K1[a], S[s].tau1[i], S[ss].tau2[j]);  // f+
						}
						if( occupied(X.state[0], ss) ){  // wind
							temp -= D0_calc_interp(K1[a], S[s].tau1[i], prm.beta);
							temp += D0_calc_interp(K1[a], S[s].tau1[i], 0);
						}
						if( a==0 )  VX_temp[s].fac_g1[i] += temp;
						else        VX_temp[s].fac_g1[i] += temp * (double)sign_spin2(s, ss) * spin_fac;
		// 				if( s == ss )  VX_temp[s].fac_g1[i] += temp;
		// 				else           VX_temp[s].fac_g1[i] -= temp;
					}
				}
			}
		}
		{
			int a = 1;  // only sp
			if( flag_g[a] ){
				for(int i=0; i<S[s].k; i++){
					double temp = 0;
					for(int j=0; j<Y.k; j++){
						temp -= D0_calc_interp(K1[a], S[s].tau1[i], Y.tau1[j]);  // S-
						temp += D0_calc_interp(K1[a], S[s].tau1[i], Y.tau2[j]);  // S+
					}
					VX_temp[s].fac_g1[i] += 2.* temp * (double)sign_spin1(s) * spin_fac;
				}
			}
		}
		
		// for self-energy (from g_z)
/*		double fac_1[S[s].k];
		for(int j=0; j<S[s].k; j++){
			fac_1[j] = 0;
			for(int ss=0; ss<2; ss++){
				double temp = 0;
				for(int i=0; i<X.k-1; i++){
					if( occupied(X.state[i+1], ss) ){
						temp -= D0_calc_interp(K1_sp, S[s].tau1[j], X.tau[i+1]);
						temp += D0_calc_interp(K1_sp, S[s].tau1[j], X.tau[i]);
					}
				}
// 				if( occupied(X.state[0], ss) ){  // wind
// 					
// 				}
				if( s == ss )  fac_1[j] += temp;
				else           fac_1[j] -= temp;
			}
		}
*/		
		for(int j=0; j<S[s].k; j++){
// 			int j_tau = tau_position(X.tau, X.k, S[s].tau1[j]);
// 			int flag_self = occupied(X.state[j_tau], (s+1)%2) ? 1 : 0;
			
			for(int i=0; i<S[s].k; i++){
// 				int j_tau = tau_position(X.tau, X.k, S[s].tau2[i]);
// 				int flag_self = occupied(X.state[j_tau], (s+1)%2) ? 1 : 0;
				
				// tau1: f-
				// tau2: f+
				double tau = S[s].tau1[j] - S[s].tau2[i];
				int i_tau;
				double g_tau = -S[s].mat_M[j][i] * (double)w_sign;
				if(tau>0){
					i_tau = (int)(tau / delta_tau);
				}
				else{
					i_tau = (int)((tau + prm.beta) / delta_tau);
					g_tau = -g_tau;
				}
				B->Gf[s][i_tau] += g_tau;
				
				if( VX_temp[s].flag_occup1[j] ){
// 					B->self[s][i_tau] += g_tau * prm.U;
					B->self[s][i_tau] += g_tau * prm_U_unrenorm;
// 					B->self[s][i_tau] += g_tau * prm_U_unrenorm * 0.5;
				}
// 				if( VX_temp[s].flag_occup2[i] ){
// 					B->self[s][i_tau] += g_tau * prm_U_unrenorm * 0.5;
// 				}
				
				B->self[s][i_tau] += g_tau * VX_temp[s].fac_g1[j];
// 				B->self[s][i_tau] += g_tau * fac_1[j] / 4.;  // sp
// 				B->self[s][i_tau] += g_tau * fac_1[j] / 4. * 3.;  // ********** isotropic
				
// 				double fac = fac_1[j];
// 				// subtract
// 				fac -= D0_calc_interp(K1_sp, S[s].tau1[j], S[s].tau1[j]);
// 				fac += D0_calc_interp(K1_sp, S[s].tau1[j], S[s].tau2[i]);
// 				B->self[s][i_tau] += g_tau * fac;
			}
		}
	}
	
}

static inline void measure_pm()
{
	//
	// chi_{+-}
	//
	double delta_tau = prm.beta / (double)N_TAU;
	double temp_m_xy = 0;  // T chi_{pm}
	
	for(int i=0; i<Y.k; i++){
		// tau1: S-
		// tau2: S+
		double tau = Y.tau2[i] - Y.tau1[i];
		int i_tau;
		
		if(tau>0){
			i_tau = (int)(tau / delta_tau);
		}
		else{
			i_tau = (int)((tau + prm.beta) / delta_tau);
		}
		double temp = 1./ D0_calc_interp(D0_pm, Y.tau1[i], Y.tau2[i]);
		B->chi_pm[i_tau] -= temp * (double)w_sign;
		temp_m_xy -= temp;
	}
	temp_m_xy *= 4./ (prm.beta * prm.beta);  // factor 4: 1/2-spin to Pauli spin
	B->m_xy_sqr += temp_m_xy * (double)w_sign;
	
	double m_rt = mz_sqr + 2.* temp_m_xy;
	B->m_rt_sqr += m_rt * (double)w_sign;
	
	
	double temp_m_xy_sqr = 0;
	for(int i=0; i<Y.k; i++){
		for(int j=i; j<Y.k; j++){
			temp_m_xy_sqr += 2./ ( D0_calc_interp(D0_pm, Y.tau1[i], Y.tau2[i]) * D0_calc_interp(D0_pm, Y.tau1[j], Y.tau2[j]) );
		}
	}
	temp_m_xy_sqr *= 16./ (prm.beta * prm.beta * prm.beta * prm.beta);  // factor 16: 1/2-spin to Pauli spin
	
	B->m_xy_pw4 += temp_m_xy_sqr * (double)w_sign;
	B->m_rt_pw4 += (mz_sqr * mz_sqr + 2.* mz_sqr * temp_m_xy + temp_m_xy_sqr )* (double)w_sign;
	
// 	for(int i=0; i<Y.k; i++){
// 		for(int j=0; j<Y.k; j++){
// 			double tau = Y.tau1[j] - Y.tau2[i];
// 			int i_tau;
// 			
// 			if(tau>0){
// 				i_tau = (int)(tau / delta_tau);
// 				B->chi_pm[i_tau] -= 1./ D0_calc_interp(D0, Y.tau1[j], Y.tau2[i]) * (double)w_sign;
// 			}
// 			else{
// 				i_tau = (int)((tau + prm.beta) / delta_tau);
// 				B->chi_pm[i_tau] -= 1./ D0_calc_interp(D0, Y.tau1[j], Y.tau2[i]) * (double)w_sign;
// 			}
// 		}
// 	}
	
}

static inline void measure_tp()
{
	if(X.k){
// 		printf(" X.tau\n");
// 		for(int i=0; i<=X.k; i++)  printf("   %d %lf\n", i, X.tau[i]);
		
		double tau1[X.k+1], tau2[2*X.k+1];
		int alpha[2*X.k];
		
		for(int i=0; i<X.k; i++){
			tau2[i] = X.tau[i];
			tau2[i+X.k] = X.tau[i] + prm.beta;
// 			alpha[i] = op2state_after[X.op[i]];
// 			alpha[i+X.k] = op2state_after[X.op[i]];
			alpha[i] = X.state[i+1];
			alpha[i+X.k] = X.state[i+1];
		}
		tau2[2*X.k] = X.tau[0] + 2.0 * prm.beta;
// 		alpha[2*X.k] = alpha[0];
		
// 		printf("measure_tp k=%d\n", X.k);
// 		for(int i=0; i<=2*X.k; i++){
// 			printf("%d %d %lf\n", i, alpha[i], tau2[i]);
// 		}
// 		printf("\n");
		
		for(int n=1; n<=N_TP; n++){  // exclude n=0
			for(int i=0; i<=X.k; i++)  tau1[i] = X.tau[i] + TP->tau[n];
// 			for(int i=0; i<=X.k; i++)  printf("   %d %lf\n", i, tau1[i]);
			
			double tau_comb[2*X.k+1];
			int alpha1[2*X.k+1], alpha2[2*X.k+1];
			
			int k=0, i1=0, i2=0;
			while( tau1[0] > tau2[i2] ) i2++;
			
			do{
				if( tau1[i1] > tau2[i2] ){
					tau_comb[k] = tau2[i2];
					alpha1[k] = alpha[i1-1];
					alpha2[k] = alpha[i2];
// 					printf(" a  %d  %d %d  %lf\n", k, alpha2[k], alpha1[k], tau_comb[k]);
					k++;
					i2++;
// 					if(i1==0){
// 						printf("error1\n");
// 						exit(0);
// 					}
				}
				else{
					tau_comb[k] = tau1[i1];
					alpha1[k] = alpha[i1];
					alpha2[k] = alpha[i2-1];
// 					printf(" b  %d  %d %d  %lf\n", k, alpha2[k], alpha1[k], tau_comb[k]);
					k++;
					i1++;
// 					if(i2==0){
// 						printf("error2\n");
// 						exit(0);
// 					}
				}
			}while( i1 <= X.k );
			
// 			printf(" n=%d k=%d   i1=%d i2=%d\n", n, k, i1, i2);
			for(int i=0; i<k-1; i++){
// 				printf("   %d  %d %d  %.3lf\n", i, alpha2[i], alpha1[i], tau_comb[i+1]);
				B->chi[alpha2[i]][alpha1[i]][n] += (tau_comb[i+1] - tau_comb[i]) * (double)w_sign;
			}
		}
	}
	else{
// 		int a = op2state[X.op[0]];
		int a = X.state[0];
		for(int i=1; i<=N_TP; i++)  B->chi[a][a][i] += prm.beta * (double)w_sign;
	}
}


static void eval_vx_u_sum()
{
	double omega_f0 = M_PI / prm.beta;
	// e^{i e2 tau_j} e^{-i e1 tau_i}
// 	for(int s=0; s<2; s++){
// 		for(int i=0; i<S[s].k; i++){
// 			// positive frequencies
// 			double x_i = S[s].tau2[i] * omega_f0;
// 			double x_j = S[s].tau1[i] * omega_f0;
// 			
// 			VX_temp[s].phase_i[N_VX1+0][i] = polar(1.0, -x_i);
// 			VX_temp[s].phase_j[N_VX1+0][i] = polar(1.0, x_j);
// 			
// 			complex<double> r_i = VX_temp[s].phase_i[N_VX1+0][i] * VX_temp[s].phase_i[N_VX1+0][i];
// 			complex<double> r_j = VX_temp[s].phase_j[N_VX1+0][i] * VX_temp[s].phase_j[N_VX1+0][i];
// 			
// 			for(int iw=1; iw<N_VX1+N_VX2; iw++){
// 				VX_temp[s].phase_i[N_VX1+iw][i] = VX_temp[s].phase_i[N_VX1+iw-1][i] * r_i;
// 				VX_temp[s].phase_j[N_VX1+iw][i] = VX_temp[s].phase_j[N_VX1+iw-1][i] * r_j;
// 			}
// 			
// 			// negative frequencies
// 			for(int iw=0; iw<N_VX1; iw++){
// 				VX_temp[s].phase_i[N_VX1-iw-1][i] = conj(VX_temp[s].phase_i[N_VX1+iw][i]);
// 				VX_temp[s].phase_j[N_VX1-iw-1][i] = conj(VX_temp[s].phase_j[N_VX1+iw][i]);
// 			}
// 		}
// 	}
	for(int s=0; s<2; s++){
		complex<double> r_i[S[s].k], r_j[S[s].k];
		for(int i=0; i<S[s].k; i++){
			VX_temp[s].phase_i[N_VX1+0][i] = polar(1.0, -S[s].tau2[i] * omega_f0);
			VX_temp[s].phase_j[N_VX1+0][i] = polar(1.0,  S[s].tau1[i] * omega_f0);
			
			r_i[i] = VX_temp[s].phase_i[N_VX1+0][i] * VX_temp[s].phase_i[N_VX1+0][i];
			r_j[i] = VX_temp[s].phase_j[N_VX1+0][i] * VX_temp[s].phase_j[N_VX1+0][i];
		}
		// positive frequencies
		for(int iw=1; iw<N_VX1+N_VX2; iw++){
			for(int i=0; i<S[s].k; i++){
				VX_temp[s].phase_i[N_VX1+iw][i] = VX_temp[s].phase_i[N_VX1+iw-1][i] * r_i[i];
				VX_temp[s].phase_j[N_VX1+iw][i] = VX_temp[s].phase_j[N_VX1+iw-1][i] * r_j[i];
			}
		}
		// negative frequencies
		for(int iw=0; iw<N_VX1; iw++){
			for(int i=0; i<S[s].k; i++){
				VX_temp[s].phase_i[N_VX1-iw-1][i] = conj(VX_temp[s].phase_i[N_VX1+iw][i]);
				VX_temp[s].phase_j[N_VX1-iw-1][i] = conj(VX_temp[s].phase_j[N_VX1+iw][i]);
			}
		}
	}
	
	// u1(e1, e2)
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
			for(int iw2=0; iw2<2*N_VX1+N_VX2; iw2++){
				VX_temp[s].u[iw1][iw2] = 0;
				VX_temp[s].u2[iw1][iw2] = 0;
				
				for(int j=0; j<S[s].k; j++){
					for(int i=0; i<S[s].k; i++){
// 						complex<double> u = S[s].mat_M[j][i] * VX_temp[s].phase_j[iw2][j] * VX_temp[s].phase_i[iw1][i];
						complex<double> u = S[s].mat_M[j][i] * VX_temp[s].phase_j[iw1][j] * VX_temp[s].phase_i[iw2][i];
						VX_temp[s].u[iw1][iw2] += u;
						
						#if IMPROV_ESTIM
						if( VX_temp[s].flag_occup1[j] ){
							VX_temp[s].u2[iw1][iw2] += u * prm_U_unrenorm;
						}
						VX_temp[s].u2[iw1][iw2] += u * VX_temp[s].fac_g1[j];
						#endif
					}
				}
				VX_temp[s].u[iw1][iw2] /= prm.beta;
				VX_temp[s].u2[iw1][iw2] /= prm.beta;
			}
		}
	}
// 	for(int s=0; s<2; s++){
// 		for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
// 			for(int iw2=0; iw2<2*N_VX1+N_VX2; iw2++){
// 				VX_temp[s].u[iw1][iw2] = VX_temp[s].u2[iw1][iw2] = 0;
// 			}
// 		}
// 		for(int j=0; j<S[s].k; j++){
// 			for(int i=0; i<S[s].k; i++){
// 				for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
// 					complex<double> u_temp = S[s].mat_M[j][i] * VX_temp[s].phase_i[iw1][i];
// 					for(int iw2=0; iw2<2*N_VX1+N_VX2; iw2++){
// 						complex<double> u = u_temp * VX_temp[s].phase_j[iw2][j];
// 						VX_temp[s].u[iw1][iw2] += u;
// 						if( VX_temp[s].flag_occup1[j] ){
// 							VX_temp[s].u2[iw1][iw2] += u * prm_U_unrenorm;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
}

// projection phi(t) : triangle, width 2a, hight 1/a
static double proj_phi_t(double tau, int nw)
{
	double a = prm.beta / double(nw);
	return ( 1. -  fabs(tau) / a ) / a;
	
	// slower than above
// 	double a_inv = double(nw) / prm.beta;
// 	return ( 1. -  fabs(tau) * a_inv ) * a_inv;
}
static double proj_phi_w(double wa)
{
	// wa = w*a
	return wa<1e-3 ? 1.- wa*wa/12. : 2./ (wa*wa) * (1.-cos(wa));
	// Taylor expansion if wa<1e-3
}

// compute l_j and phi_j using tau_j
//  tau_j[K]     : set of tau of the data M(tau, tau')
//  l_j[K][2M]   : poisition onto which M(tau) is projected
//  phi_j[K][2M] : weight
static void vx_proj(double *tau_j, int **l_j, double **phi_j, int K, int nw, int M)
{
	int N = nw * M;
	double delta_tau = prm.beta / double(N);
	
	for(int i=0; i<K; i++){
		double tau = tau_j[i];
		int p_tau = (int)(tau / delta_tau);  // position of tau
		
		for(int m=0; m<2*M; m++){
			int l = p_tau + M - m;  // (p-M:p+M]
			l_j[i][m] = (l+N)%N;  // [0:N)
			
// 			double dist = fabs(tau / delta_tau - double(l));  // [0:M)
// 			phi_j[i][m] = proj_phi_t(dist/double(M), nw);
			
			double dist = fabs(tau - double(l) * delta_tau);
			phi_j[i][m] = proj_phi_t(dist, nw);
			
			// impose anti-periodicity for M
			if( l<0 || l>=N )  phi_j[i][m] = -phi_j[i][m];
		}
	}
}
/*
static void fft_M(double **M_t, complex<double> **M_w, int N)
{
	static fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N*N);
	static fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N*N);
	static fftw_plan p = fftw_plan_dft_2d(2*N, 2*N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	// FFTW_BACKWARD : exp(+iwt)
	
	// extend period from beta to 2*beta
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			{
				int k = i*(2*N) + j;
				in[k][0] = M_t[i][j];
				in[k][1] = 0;
			}
			{
				int k = (i+N)*(2*N) + j;
				in[k][0] = -M_t[i][j];
				in[k][1] = 0;
			}
			{
				int k = i*(2*N) + (j+N);
				in[k][0] = -M_t[i][j];
				in[k][1] = 0;
			}
			{
				int k = (i+N)*(2*N) + (j+N);
				in[k][0] = M_t[i][j];
				in[k][1] = 0;
			}
		}
	}
	
	fftw_execute(p);
	
	double fac = prm.beta / double(4*N*N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			int mj = N - 1 - j;  // w(j) -> -w(j)
			int k = (2*i+1)*(2*N) + (2*mj+1);  // pick up odd (fermionic) frequencies
			M_w[i][j] = complex<double>( out[k][0], out[k][1] ) * fac;
		}
	}
	
// 	fftw_destroy_plan(p);
// 	fftw_free(in);
// 	fftw_free(out);
}
*/
static void fft_M(double **M_t, complex<double> **M_wp, complex<double> **M_wm, int N, int NW_p)
{
	static complex<double> *temp = (complex<double>*) fftw_malloc(sizeof(complex<double>) * N*NW_p);
	
	//
	// M(j,i) -> M(j,w')
	// Real FFT:
	// N real numers to N/2 complex numbers for each j (j runs N values)
	//
	{
		static double *in = (double*) fftw_malloc(sizeof(double) * 2*N);
		static fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+1));
		static fftw_plan p1 = fftw_plan_dft_r2c_1d(2*N, in, out, FFTW_ESTIMATE);
		// FFTW_FORWARD : exp(-iwt)
		
		for(int j=0; j<N; j++){
			// extend period from beta to 2*beta
			for(int i=0; i<N; i++){
				in[i] = M_t[j][i];
				in[i+N] = -M_t[j][i];
			}
			
			fftw_execute(p1);
			
			for(int i=0; i<NW_p; i++){  // copy only NW_p
				int w = 2*i+1;  // pick up odd (fermionic) frequencies
				temp[i*N+j] = complex<double>( out[w][0], out[w][1] );  // [i][j]
			}
		}
// 		fftw_destroy_plan(p);
// 		fftw_free(in);
// 		fftw_free(out);
	}
	
	//
	// M(j,w') -> M(w,w')
	// Complex FFT:
	// N complex numers to N complex numbers for each w' (w' runs NW_p << N values)
	//
	{
		static fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
		static fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
		static fftw_plan p2 = fftw_plan_dft_1d(2*N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		// FFTW_BACKWARD : exp(+iwt)
		
		double fac = prm.beta / double(4*N*N);
		
		for(int i=0; i<NW_p; i++){  // only NW_p (save time)
			for(int j=0; j<N; j++){
				in[j][0] = real(temp[i*N+j]);  // [i][j]
				in[j][1] = imag(temp[i*N+j]);
				in[j+N][0] = -real(temp[i*N+j]);
				in[j+N][1] = -imag(temp[i*N+j]);
			}
			
			fftw_execute(p2);
			
			for(int j=0; j<NW_p; j++){
				// pick up odd (fermionic) frequencies
				int wp = 2*j+1;  // w>0
				int wm = 2*(N-1-j)+1;  // w<0
				M_wp[j][i] = complex<double>( out[wp][0], out[wp][1] ) * fac;
				M_wm[j][i] = complex<double>( out[wm][0], out[wm][1] ) * fac;
			}
		}
	}
}
static void eval_vx_u_fft()
{
// 	const int NW = 2*N_VX1+N_VX2;  // # of freq required
	const int NW = 2*(N_VX1+N_VX2);  // # of freq required
	const int NW_p = N_VX1+N_VX2;    //  - positive freq
	const int NW_m = N_VX1;          //  - negative freq
	const int N = NW * M_OS;         // # of freq computed
	
	const int A = IMPROV_ESTIM ? 2 : 1;
	
	static double ***M_t;  // [A][N][N],  M_t[1] for improved estimator
	static complex<double> **M_wp;  // (w>0, w'>0)
	static complex<double> **M_wm;  // (w<0, w'>0)
	static int **l_j, **l_i;
	static double **phi_j, **phi_i;
	
	static double fac[NW_p][NW_p];
	
	static bool flag_first_time = true;
	if( flag_first_time ){  // execute only once
		// allocate memory
		array3_alloc(M_t, A, N, N);
// 		array2_alloc(M_w, N, N);
		array2_alloc(M_wp, NW_p, NW_p);
		array2_alloc(M_wm, NW_p, NW_p);
		array2_alloc(l_j, N_K, 2*M_OS);
		array2_alloc(l_i, N_K, 2*M_OS);
		array2_alloc(phi_j, N_K, 2*M_OS);
		array2_alloc(phi_i, N_K, 2*M_OS);
		
		// renormalize factor
		double phi_w[NW_p];
		for(int i=0; i<NW_p; i++){
			double w = double(2*i+1) / double(NW) * M_PI;
			phi_w[i] = proj_phi_w(w);
		}
		for(int i=0; i<NW_p; i++){
			for(int j=0; j<NW_p; j++)  fac[i][j] = 1./ (phi_w[i] * phi_w[j]);
		}
		flag_first_time = false;
	}
	
	for(int s=0; s<2; s++){
		array3_zero(M_t, A, N, N);
		
		// prepare projection factor phi and list of indices l
		// to create linear mesh from tau1[] and tau2[]
		vx_proj(S[s].tau1, l_j, phi_j, S[s].k, NW, M_OS);
		vx_proj(S[s].tau2, l_i, phi_i, S[s].k, NW, M_OS);
		
		// project M onto linear mesh
		for(int j=0; j<S[s].k; j++){
			for(int i=0; i<S[s].k; i++){
				for(int m1=0; m1<2*M_OS; m1++){
					for(int m2=0; m2<2*M_OS; m2++){
						int l1 = l_j[j][m1];
						int l2 = l_i[i][m2];
// 						M_t[l1][l2] += phi_j[j][m1] * S[s].mat_M[j][i] * phi_i[i][m2];
						
						double u = phi_j[j][m1] * S[s].mat_M[j][i] * phi_i[i][m2];
						M_t[0][l1][l2] += u;
						
						#if IMPROV_ESTIM
						if( VX_temp[s].flag_occup1[j] ){
							M_t[1][l1][l2] += u * prm_U_unrenorm;
						}
						M_t[1][l1][l2] += u * VX_temp[s].fac_g1[j];
						#endif
					}
				}
			}
		}
		
		// pointer to VX_temp[s].u[2*N_VX1+N_VX2][2*N_VX1+N_VX2]
		complex<double> (*p_u[2])[2*N_VX1+N_VX2] = {VX_temp[s].u, VX_temp[s].u2};
		for(int a=0; a<A; a++){
			// a=1 for improved estimator
			
			// FFT
			fft_M(M_t[a], M_wp, M_wm, N, NW_p);
			
			// renormalize
			for(int j=0; j<NW_p; j++){
				for(int i=0; i<NW_p; i++){
					M_wp[j][i] *= fac[j][i];
					M_wm[j][i] *= fac[j][i];
				}
			}
			
			// copy
			for(int j=0; j<NW_p; j++){  // w1>0
				for(int i=0; i<NW_p; i++){  // w2>0
					p_u[a][N_VX1+j][N_VX1+i] = M_wp[j][i];
				}
				for(int i=0; i<NW_m; i++){  // w2<0
					// M(w1>0, w2<0) = M*(-w1<0, -w2>0)
					p_u[a][N_VX1+j][N_VX1-1-i] = conj(M_wm[j][i]);
				}
			}
			for(int j=0; j<NW_m; j++){  // w1<0
				for(int i=0; i<NW_p; i++){  // w2>0
					// M(w1<0, w2>0)
					p_u[a][N_VX1-1-j][N_VX1+i] = M_wm[j][i];
				}
				for(int i=0; i<NW_m; i++){  // w2<0
					// M(w1<0, w2<0) = M*(-w1>0, -w2>0)
					p_u[a][N_VX1-1-j][N_VX1-1-i] = conj(M_wp[j][i]);
				}
			}
		}
	}
	
	// free memory
// 	array2_free(M_t, N);
// 	array2_free(M_w, N);
// 	array2_free(l_j, N_K);
// 	array2_free(l_i, N_K);
// 	array2_free(phi_j, N_K);
// 	array2_free(phi_i, N_K);
}

// for test
static void compare_vx_u()
{
	eval_vx_u_sum();
	
	// copy
	complex<double> ***u_fft;
	complex<double> ***u2_fft;
	array3_alloc(u_fft, 2, 2*N_VX1+N_VX2, 2*N_VX1+N_VX2);
	array3_alloc(u2_fft, 2, 2*N_VX1+N_VX2, 2*N_VX1+N_VX2);
	{
		for(int s=0; s<2; s++){
			for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
				for(int iw2=0; iw2<2*N_VX1+N_VX2; iw2++){
					u_fft[s][iw1][iw2] = VX_temp[s].u[iw1][iw2];
					u2_fft[s][iw1][iw2] = VX_temp[s].u2[iw1][iw2];
				}
			}
		}
	}
	
	eval_vx_u_fft();
	
	// check accuracy of N-FFT
	{
		double max_err[2] = {0}, ave_err[2] = {0};
		for(int s=0; s<2; s++){
			for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
				for(int iw2=0; iw2<2*N_VX1+N_VX2; iw2++){
					double err[2];
					err[0] = abs( VX_temp[s].u[iw1][iw2] - u_fft[s][iw1][iw2] );
					err[1] = abs( VX_temp[s].u2[iw1][iw2] - u2_fft[s][iw1][iw2] );
					for(int a=0; a<2; a++){
						max_err[a] = max(max_err[a], err[a]);
						ave_err[a] += err[a];
					}
				}
			}
		}
		printf("\n accuracy of N-FFT:\n");
		for(int a=0; a<2; a++){
			ave_err[a] /= pow(2*N_VX1+N_VX2, 2);
			printf("  error of u%d:  max = %.2e  ave = %.2e\n", a, max_err[a], ave_err[a]);
		}
		
		FILE *fp=fopen("test_vx.dat", "w");
// 		int iw2=0;
		int iw2=N_VX1;
		for(int iw1=0; iw1<2*N_VX1+N_VX2; iw1++){
			fprintf(fp, "%d %d", iw1-N_VX1, iw2-N_VX1);
			fprintf(fp, " %.4lf %.4lf", real(u_fft[0][iw1][iw2]), imag(u_fft[0][iw1][iw2]));
			fprintf(fp, " %.4lf %.4lf", real(VX_temp[0].u[iw1][iw2]), imag(VX_temp[0].u[iw1][iw2]));
			fprintf(fp, " %.4lf %.4lf", real(u2_fft[0][iw1][iw2]), imag(u2_fft[0][iw1][iw2]));
			fprintf(fp, " %.4lf %.4lf", real(VX_temp[0].u2[iw1][iw2]), imag(VX_temp[0].u2[iw1][iw2]));
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" 'test_vx.dat\n");
	}
	
	array3_free(u_fft, 2, 2*N_VX1+N_VX2);
	array3_free(u2_fft, 2, 2*N_VX1+N_VX2);
}
static inline void measure_vx()
{
	// compute u and u2
	//   VX_temp[s].u[2*N_VX1+N_VX2][2*N_VX1+N_VX2]
	//   VX_temp[s].u2[2*N_VX1+N_VX2][2*N_VX1+N_VX2]
	//
	//  freq = ( 2 * (iw-N_VX1) + 1 ) * pi * T
	//  (N_VX1 negative freq, N_VX1+N_VX2 positive freq )
	
	#if VX_NFFT == 0  // direct sum (no FFT)
		eval_vx_u_sum();
	
	#elif VX_NFFT == 1  // Non-equidistant FFT
		eval_vx_u_fft();
	
	#elif VX_NFFT == 2  // TEST: compare the above two
		compare_vx_u();
		exit(0);
		
	#endif  // VX_NFFT
	
	
	// two-particle Green function
	// charge and longtudinal spin
	for(int s1=0; s1<2; s1++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G4_lo[s1][s1][iw1][iw2][nu] -= VX_temp[s1].u[iw1][iw2] * VX_temp[s1].u[iw2+nu][iw1+nu] * (double)w_sign;
					#if IMPROV_ESTIM
// 					B->vx_gamma_lo[s1][s1][iw1][iw2][nu] -= VX_temp[s1].u2[iw1][iw2] * VX_temp[s1].u[iw2+nu][iw1+nu] * (double)w_sign;
					B->vx_gamma_lo[s1][s1][iw1][iw2][nu] -= VX_temp[s1].u[iw1][iw2] * VX_temp[s1].u2[iw2+nu][iw1+nu] * (double)w_sign;
					#endif
				}
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						B->vx_G4_lo[s1][s2][iw1][iw2][nu] += VX_temp[s1].u[iw1][iw1+nu] * VX_temp[s2].u[iw2+nu][iw2] * (double)w_sign;
						#if IMPROV_ESTIM
// 						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] += VX_temp[s1].u2[iw1][iw1+nu] * VX_temp[s2].u[iw2+nu][iw2] * (double)w_sign;
						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] += VX_temp[s1].u[iw1][iw1+nu] * VX_temp[s2].u2[iw2+nu][iw2] * (double)w_sign;
						#endif
					}
				}
			}
		}
	}
	// transverse spin
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G4_tr[s][iw1][iw2][nu] -= VX_temp[s].u[iw1][iw2] * VX_temp[(s+1)%2].u[iw2+nu][iw1+nu] * (double)w_sign;
					#if IMPROV_ESTIM
// 					B->vx_gamma_tr[s][iw1][iw2][nu] -= VX_temp[s].u2[iw1][iw2] * VX_temp[(s+1)%2].u[iw2+nu][iw1+nu] * (double)w_sign;
					B->vx_gamma_tr[s][iw1][iw2][nu] -= VX_temp[s].u[iw1][iw2] * VX_temp[(s+1)%2].u2[iw2+nu][iw1+nu] * (double)w_sign;
					#endif
				}
			}
		}
	}
	
	// single-particle Green function
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			B->vx_Gf[s][iw] -= VX_temp[s].u[iw][iw] * (double)w_sign;
			B->vx_self[s][iw] -= VX_temp[s].u2[iw][iw] * (double)w_sign;
		}
	}
	
	//
	// Three-point vertex
	//
	double omega_b = 2. * M_PI / prm.beta;
	
	complex<double> temp_init[N_VX2] = {0};
	for(int i=0; i<Y.k; i++){
		complex<double> phase_1 = polar(1.0,  Y.tau1[i] * omega_b);
		complex<double> phase_2 = polar(1.0,  Y.tau2[i] * omega_b);
		complex<double> phase_1_nu = 1.;
		complex<double> phase_2_nu = 1.;
		for(int nu=1; nu<N_VX2; nu++){
			// tau1: S-
			// tau2: S+
			phase_1_nu *= phase_1;
			phase_2_nu *= phase_2;
			temp_init[nu] += phase_1_nu;
			temp_init[nu] -= phase_2_nu;
		}
	}
	for(int s=0; s<2; s++){
		complex<double> temp[N_VX2];
		for(int nu=1; nu<N_VX2; nu++)  temp[nu] = temp_init[nu] * (double)sign_spin1(s);
		// nu=0 is computed in function measure_sp
		
		for(int i=0; i<S[s].k; i++){
			complex<double> phase_1 = polar(1.0,  S[s].tau1[i] * omega_b);
			complex<double> phase_2 = polar(1.0,  S[s].tau2[i] * omega_b);
			complex<double> phase_1_nu = 1.;
			complex<double> phase_2_nu = 1.;
			for(int nu=1; nu<N_VX2; nu++){
				// tau1: f-
				// tau2: f+
				phase_1_nu *= phase_1;
				phase_2_nu *= phase_2;
				temp[nu] += phase_1_nu;
				temp[nu] -= phase_2_nu;
			}
		}
		for(int nu=1; nu<N_VX2; nu++){
			VX_temp[s].occup[nu] = temp[nu] / (IMAG * double(nu) * omega_b);
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G3[s1][s2][iw][nu] -= VX_temp[s1].u[iw][iw+nu] * VX_temp[s2].occup[nu] * (double)w_sign;
				}
			}
		}
	}
	// ********************FFT
	
	// susceptibility
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int nu=0; nu<N_VX2; nu++){
				B->vx_suscep[s1][s2][nu] += VX_temp[s1].occup[nu] * conj(VX_temp[s2].occup[nu]) * (double)w_sign;
			}
		}
	}
}

static void func_measure0()
{
	measure_stat();
}

static void func_measure1()
{
	measure_stat();
	measure_sp();
	measure_pm();
}

static void func_measure2()
{
	measure_stat();
	measure_sp();
	measure_pm();
	measure_tp();
	
// 	#if CHI_TR
// // 	measure_tp_tr1();
// 	measure_tp_tr2();
// 	#endif // CHI_TR
}

static void func_measure3()
{
	measure_stat();
	measure_sp();
	measure_pm();
	measure_tp();
	measure_vx();
}



//============================================================================
//
// AVERAGE DATA IN A BIN
//

static void averagebin_stat(int n_sample)
{
	double B_ave_sign = (double)B->ave_sign / (double)n_sample;
	PQ.ave_sign += B_ave_sign;
	PQ.ave_sign_err += pow(B_ave_sign, 2);
	
	double sum_n_k = 0;
	for(int i=0; i<N_K; i++)  sum_n_k += (double)B->n_k[0][i];
	
	double B_ave_Z_V[2] = {0}, B_ave_Z_g = 0;
	double B_n_l[N_K];
	for(int i=0; i<N_K; i++){
		for(int s=0; s<2; s++){
			double B_n_k = (double)B->n_k[s][i] / sum_n_k;
			PQ.Z_V[s][i] += B_n_k;
			PQ.Z_V_err[s][i] += pow(B_n_k, 2);
			
			B_ave_Z_V[s] += B_n_k * (double)i;
		}
		
		B_n_l[i] = (double)B->n_l[i] / sum_n_k;
		PQ.Z_g[i] += B_n_l[i];
		PQ.Z_g_err[i] += pow(B_n_l[i], 2);
		
		B_ave_Z_g += B_n_l[i] * (double)i;
		
	}
	
	for(int s=0; s<2; s++){
		PQ.ave_Z_V[s] += B_ave_Z_V[s];
		PQ.ave_Z_V_err[s] += pow(B_ave_Z_V[s], 2);
	}
	PQ.ave_Z_g += B_ave_Z_g;
	PQ.ave_Z_g_err += pow(B_ave_Z_g, 2);
	
	
	double temp_free_energy = 0;
	for(int i=0; i<N_K; i++){
		temp_free_energy += B_n_l[i] * exp( WL.log_gk[i] - WL.log_gk[0] );
	}
	double B_free_energy = ( log(B_n_l[0]) - log(temp_free_energy) ) / prm.beta;
	PQ.free_energy += B_free_energy;
	PQ.free_energy_err += pow(B_free_energy, 2);
}

static void averagebin_sp(int n_sample)
{
	for(int a=0; a<4; a++){
		B->occup[a] /= prm.beta * (double)n_sample;
	}
	{
		double B_occup_ud[2] = {
			B->occup[UP] + B->occup[DBL], 
			B->occup[DW] + B->occup[DBL]
		};
		for(int s=0; s<2; s++){
			PQ.occup_ud[s] += B_occup_ud[s];
			PQ.occup_ud_err[s] += pow(B_occup_ud[s], 2);
		}
		
		double occup_tot = B_occup_ud[0] + B_occup_ud[1];
		PQ.occup_tot += occup_tot;
		PQ.occup_tot_err += pow(occup_tot, 2);
		
		double occup_dif = B_occup_ud[0] - B_occup_ud[1];
		PQ.occup_dif += occup_dif;
		PQ.occup_dif_err += pow(occup_dif, 2);
		
		PQ.occup_dbl += B->occup[DBL];
		PQ.occup_dbl_err += pow(B->occup[DBL], 2);
	}
	
	B->m_z_sqr /= (double)n_sample;
	PQ.m_z_sqr += B->m_z_sqr;
	PQ.m_z_sqr_err += pow(B->m_z_sqr, 2);
	
	B->m_z_pw4 /= (double)n_sample;
	PQ.m_z_pw4 += B->m_z_pw4;
	PQ.m_z_pw4_err += pow(B->m_z_pw4, 2);
	
	double delta_tau = prm.beta / (double)N_TAU;
	for(int s=0; s<2; s++){
		for(int i=0; i<N_TAU; i++){
			B->Gf[s][i] /= prm.beta * delta_tau * (double)n_sample;
			B->self[s][i] /= prm.beta * delta_tau * (double)n_sample;
		}
		
		for(int i=N_TAU-1; i>0; i--){
			B->Gf[s][i] = (B->Gf[s][i] + B->Gf[s][i-1]) * 0.5;
			B->self[s][i] = (B->self[s][i] + B->self[s][i-1]) * 0.5;
		}
		
// 		B->Gf[s][0] = B->f_number[s] - 1.0;
// 		B->Gf[s][N_TAU] = - B->f_number[s];
		
// 		for(int i=0; i<=N_TAU; i++){
		for(int i=1; i<N_TAU; i++){
			SP[s].Gf_tau[i] += B->Gf[s][i];
			SP[s].Gf_tau_err[i] += pow(B->Gf[s][i], 2);
			
			SP[s].G_self[i] += B->self[s][i];
// 			SP[s].G_self[i] += B->self[s][i] / 2.;
		}
	}
	
}

static void averagebin_pm(int n_sample)
{
	//
	// chi_{+-}
	//
	double delta_tau = prm.beta / (double)N_TAU;
	
	double B_stat_suscep_pm = 0;
	for(int i=0; i<N_TAU; i++){
		B->chi_pm[i] /= prm.beta * delta_tau * (double)n_sample;
		B->chi_pm[i] *= 4.* fac_chi[2];  // 1/2-spin to Pauli
		B_stat_suscep_pm += B->chi_pm[i];
	}
	for(int i=N_TAU-1; i>0; i--){
		B->chi_pm[i] = (B->chi_pm[i] + B->chi_pm[i-1]) * 0.5;
	}
	for(int i=1; i<N_TAU; i++){
		TP->chi_pm_tau[i] += B->chi_pm[i];
		TP->chi_pm_tau_err[i] += pow(B->chi_pm[i], 2);
	}
// 	double chi_tau0 = B->occup[UP] / 2.;
// 	double chi_tau1 = B->occup[DW] / 2.;
	double chi_tau0 = B->occup[UP] * 2. * fac_chi[2];
	double chi_tau1 = B->occup[DW] * 2. * fac_chi[2];
	TP->chi_pm_tau[0] += chi_tau0;
	TP->chi_pm_tau_err[0] += pow(chi_tau0, 2);
	TP->chi_pm_tau[N_TAU] += chi_tau1;
	TP->chi_pm_tau_err[N_TAU] += pow(chi_tau1, 2);
	
	// static
	B_stat_suscep_pm *= delta_tau;
	PQ.stat_suscep_pm += B_stat_suscep_pm;
	PQ.stat_suscep_pm_err += pow(B_stat_suscep_pm, 2);
	
	B->m_xy_sqr /= (double)n_sample;
	PQ.m_xy_sqr += B->m_xy_sqr;
	PQ.m_xy_sqr_err += pow(B->m_xy_sqr, 2);
	
	B->m_xy_pw4 /= (double)n_sample;
	PQ.m_xy_pw4 += B->m_xy_pw4;
	PQ.m_xy_pw4_err += pow(B->m_xy_pw4, 2);
	
	B->m_rt_sqr /= (double)n_sample;
	PQ.m_rt_sqr += B->m_rt_sqr;
	PQ.m_rt_sqr_err += pow(B->m_rt_sqr, 2);
	
	B->m_rt_pw4 /= (double)n_sample;
	PQ.m_rt_pw4 += B->m_rt_pw4;
	PQ.m_rt_pw4_err += pow(B->m_rt_pw4, 2);
}

static void averagebin_tp(int n_sample)
{
	for(int a1=0; a1<4; a1++){
		for(int a2=0; a2<4; a2++){
			for(int i=1; i<=N_TP; i++){
				B->chi[a1][a2][i] /= prm.beta * (double)n_sample;
			}
			
			if(a1==a2)  B->chi[a1][a2][0] = B->occup[a1];
			else        B->chi[a1][a2][0] = 0;
			
			double temp = B->occup[a1] * B->occup[a2];
			for(int i=0; i<=N_TP; i++)  B->chi[a1][a2][i] -= temp;
		}
	}
	
	double B_chi_sp[N_TP+1], B_chi_ch[N_TP+1];
	for(int i=0; i<=N_TP; i++){
// 		double diag = B->chi[UP][UP][i] + B->chi[DW][DW][i];
// 		double offd = B->chi[UP][DW][i] + B->chi[DW][UP][i];
		
		double B_chi_uu = B->chi[UP][UP][i] + B->chi[UP][DBL][i] + B->chi[DBL][UP][i] + B->chi[DBL][DBL][i];
		double B_chi_ud = B->chi[UP][DW][i] + B->chi[UP][DBL][i] + B->chi[DBL][DW][i] + B->chi[DBL][DBL][i];
		double B_chi_du = B->chi[DW][UP][i] + B->chi[DW][DBL][i] + B->chi[DBL][UP][i] + B->chi[DBL][DBL][i];
		double B_chi_dd = B->chi[DW][DW][i] + B->chi[DW][DBL][i] + B->chi[DBL][DW][i] + B->chi[DBL][DBL][i];
		double diag = B_chi_uu + B_chi_dd;
		double offd = B_chi_ud + B_chi_du;
		
// 		B_chi_sp[i] = (diag - offd) / 4.;  // < S_z ; S_z >
		B_chi_sp[i] = (diag - offd) * fac_chi[1];
		B_chi_ch[i] = (diag + offd) * fac_chi[0];
		
		TP->chi_sp_tau[i] += B_chi_sp[i];
		TP->chi_ch_tau[i] += B_chi_ch[i];
		TP->chi_sp_tau_err[i] += pow(B_chi_sp[i], 2);
		TP->chi_ch_tau_err[i] += pow(B_chi_ch[i], 2);
	}
	
	// static susceptibility
	//
	// trapezoidal rule (non-linear mesh)
	double B_stat_suscep_sp = 0;
	double B_stat_suscep_ch = 0;
	for(int i=0; i<N_TP; i++){
		B_stat_suscep_sp += (B_chi_sp[i+1] + B_chi_sp[i]) * (TP->tau[i+1] - TP->tau[i]);
		B_stat_suscep_ch += (B_chi_ch[i+1] + B_chi_ch[i]) * (TP->tau[i+1] - TP->tau[i]);
	}
	// extend to [0:beta] (factor 0.5 was omitted instead)
	
	PQ.stat_suscep_sp += B_stat_suscep_sp;
	PQ.stat_suscep_ch += B_stat_suscep_ch;
	PQ.stat_suscep_sp_err += pow(B_stat_suscep_sp, 2);
	PQ.stat_suscep_ch_err += pow(B_stat_suscep_ch, 2);
}

static void averagebin_vx(int n_sample)
{
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						B->vx_G4_lo[s1][s2][iw1][iw2][nu] /= (double)n_sample;
						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] /= (double)n_sample;
						VX->Gfour_lo[s1][s2][iw1][iw2][nu] += B->vx_G4_lo[s1][s2][iw1][iw2][nu];
						VX->gamma_lo[s1][s2][iw1][iw2][nu] += B->vx_gamma_lo[s1][s2][iw1][iw2][nu];
					}
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G4_tr[s][iw1][iw2][nu] /= (double)n_sample;
					B->vx_gamma_tr[s][iw1][iw2][nu] /= (double)n_sample;
					VX->Gfour_tr[s][iw1][iw2][nu] += B->vx_G4_tr[s][iw1][iw2][nu];
					VX->gamma_tr[s][iw1][iw2][nu] += B->vx_gamma_tr[s][iw1][iw2][nu];
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			B->vx_Gf[s][iw] /= (double)n_sample;
			B->vx_self[s][iw] /= (double)n_sample;
			VX->G[s][iw] += B->vx_Gf[s][iw];
			VX->self[s][iw] += B->vx_self[s][iw];
		}
	}
	
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int nu=0; nu<N_VX2; nu++){
				B->vx_suscep[s1][s2][nu] /= (double)n_sample;
				VX->suscep_lo[s1][s2][nu] += B->vx_suscep[s1][s2][nu] / prm.beta;
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					B->vx_G3[s1][s2][iw][nu] /= (double)n_sample;
					VX->Gthree_lo[s1][s2][iw][nu] += B->vx_G3[s1][s2][iw][nu];
				}
			}
		}
	}
}
/*
static void averagebin_vx(int n_sample)
{
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
// 						B->vx_G4_lo[s1][s2][iw1][iw2][nu] /= prm.beta * prm.beta * (double)n_sample;
// 						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] /= prm.beta * prm.beta * (double)n_sample;
						B->vx_G4_lo[s1][s2][iw1][iw2][nu] /= prm.beta * (double)n_sample;  // ******
						B->vx_gamma_lo[s1][s2][iw1][iw2][nu] /= prm.beta * (double)n_sample;
					}
				}
			}
		}
	}
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				complex<double> diag, offd;
				diag = B->vx_G4_lo[0][0][iw1][iw2][nu] + B->vx_G4_lo[1][1][iw1][iw2][nu];
				offd = B->vx_G4_lo[0][1][iw1][iw2][nu] + B->vx_G4_lo[1][0][iw1][iw2][nu];
				VX->chi_ch[iw1][iw2][nu] += (diag + offd) / 2.;
				VX->chi_sp[iw1][iw2][nu] += (diag - offd) / 2.;
				
				diag = B->vx_gamma_lo[0][0][iw1][iw2][nu] + B->vx_gamma_lo[1][1][iw1][iw2][nu];
				offd = B->vx_gamma_lo[0][1][iw1][iw2][nu] + B->vx_gamma_lo[1][0][iw1][iw2][nu];
				VX->vx_ch[iw1][iw2][nu] += (diag + offd) / 2.;
				VX->vx_sp[iw1][iw2][nu] += (diag - offd) / 2.;
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			B->vx_Gf[s][iw] /= prm.beta * (double)n_sample;
			B->vx_self[s][iw] /= prm.beta * (double)n_sample;
			VX->G[s][iw] += B->vx_Gf[s][iw];
			VX->self[s][iw] += B->vx_self[s][iw];
		}
	}
}
*/

static void func_averagebin0(int n_sample)
{
	averagebin_stat(n_sample);
}
static void func_averagebin1(int n_sample)
{
	averagebin_stat(n_sample);
	averagebin_sp(n_sample);
	averagebin_pm(n_sample);
}
static void func_averagebin2(int n_sample)
{
	averagebin_stat(n_sample);
	averagebin_sp(n_sample);
	averagebin_pm(n_sample);
	averagebin_tp(n_sample);
}
static void func_averagebin3(int n_sample)
{
	averagebin_stat(n_sample);
	averagebin_sp(n_sample);
	averagebin_pm(n_sample);
	averagebin_tp(n_sample);
	averagebin_vx(n_sample);
}


//============================================================================
//
// AVERAGE DATA TO ESTIMATE ERROR
//

static inline void average_sub(double &msr, double &msr_err, int n_bin, double sign)
{
	msr /= (double)n_bin;
	msr_err /= (double)n_bin;
	msr_err = sqrt(msr_err - msr*msr) * sqrt( (double)n_bin / (double)(n_bin-1) );
	
	msr /= sign;
	msr_err /= sign;
}
static inline void average_sub(complex<double> &msr, complex<double> &msr_err, int n_bin, double sign)
{
	msr /= (double)n_bin;
	msr_err /= (double)n_bin;
	
	double fac_err = sqrt( (double)n_bin / (double)(n_bin-1) );
	msr_err = sqrt(real(msr_err) - pow(real(msr), 2)) * fac_err
	        + sqrt(imag(msr_err) - pow(imag(msr), 2)) * fac_err * IMAG;
	
	msr /= sign;
	msr_err /= sign;
}

static inline void average_sub(double &msr, int n_bin, double sign)
{
	double dummy = 0;
	average_sub(msr, dummy, n_bin, sign);
}
static inline void average_sub(complex<double> &msr, int n_bin, double sign)
{
	complex<double> dummy = 0;
	average_sub(msr, dummy, n_bin, sign);
}

static void fft_chi_after_interp(double *chi_tau, complex<double> *chi_omega)
{
	double chi_tau2[2*N_TP2+1];
	
	// interpolation
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (INTERP, N_TP+1);
	gsl_spline_init(spline, TP->tau, chi_tau, N_TP+1);
	
	for(int i=0; i<=N_TP2; i++){
		double TP2_tau = (double)i * prm.beta / (double)(2*N_TP2);
		chi_tau2[i] = gsl_spline_eval(spline, TP2_tau, acc);
	}
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	
	// extend range from [0:beta/2] to [0:beta]
	for(int i=0; i<N_TP2; i++)  chi_tau2[2*N_TP2-i] = chi_tau2[i];
	fft_boson_radix2_tau2omega(chi_tau2, chi_omega, prm.beta, 2*N_TP2);
}

/*
void analys_stat(unsigned long *n_k, double *d_n_k, int N, char *str)
{
	double sum_n_k = 0;
	for(int i=0; i<N; i++)  sum_n_k += (double)n_k[i];
	
	for(int i=0; i<N; i++)  d_n_k[i] = (double)n_k[i] / sum_n_k;
	
	double ave_k = 0;
	for(int i=0; i<N; i++){
		ave_k += (double)i * d_n_k[i];
	}
	
	int peak_nk = 0;
	for(int i=1; i<N; i++){
		if( n_k[peak_nk] < n_k[i] )  peak_nk = i;
	}
	
	int max_nk = 0;
	for(int i=1; i<N; i++){
		if( n_k[i] != 0 )  max_nk = i;
	}
	
	if(DISPLAY){
		printf("%s", str);
		printf(" ave %6.2lf", ave_k);
		printf(" | peak %3d", peak_nk);
		printf(" | max %3d\n", max_nk);
	}
	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "%s", str);
	fprintf(fp_log, " ave %6.2lf", ave_k);
	fprintf(fp_log, " | peak %3d", peak_nk);
	fprintf(fp_log, " | max %3d\n", max_nk);
	fclose(fp_log);
}
*/
/*
void analys_stat(double *d_n_k, double ave_k, int N, char *str)
{
	int max_nk = 0;
	for(int i=1; i<N; i++){
		if( d_n_k[i] > 0 )  max_nk = i;
	}
	
	if(DISPLAY){
		printf("%s", str);
		printf(" ave %6.8lf | max %3d\n", ave_k, max_nk);
	}
	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "%s", str);
	fprintf(fp_log, " ave %6.8lf | max %3d\n", ave_k, max_nk);
	fclose(fp_log);
}
*/

static int max_nonzero(double *d_n_k, int N)
{
	int max_nk = 0;
	for(int i=1; i<N; i++){
		if( d_n_k[i] > 0 )  max_nk = i;
	}
	return max_nk;
}

static void average_stat(int n_bin)
{
	average_sub(PQ.ave_sign, PQ.ave_sign_err, n_bin, 1.0);
	
	if(DISPLAY){
		printf("\n Average sign\n");
		printf("   %.6lf +- %.6lf", PQ.ave_sign, PQ.ave_sign_err);
		printf("  (deviation from 1:");
		printf(" %.3e +- %.3e)\n", 1.0-PQ.ave_sign, PQ.ave_sign_err);
	}
	
	for(int i=0; i<N_K; i++){
		for(int s=0; s<2; s++){
			average_sub(PQ.Z_V[s][i], PQ.Z_V_err[s][i], n_bin, PQ.ave_sign);
		}
		average_sub(PQ.Z_g[i], PQ.Z_g_err[i], n_bin, PQ.ave_sign);
	}
	
	for(int s=0; s<2; s++){
		average_sub(PQ.ave_Z_V[s], PQ.ave_Z_V_err[s], n_bin, PQ.ave_sign);
	}
	average_sub(PQ.ave_Z_g, PQ.ave_Z_g_err, n_bin, PQ.ave_sign);
	
	if(DISPLAY){
		printf("\n Expansion order\n");
		for(int s=0; s<2; s++){
			printf("  V^2k %d  : %6.2lf +- %.2lf  (max %d)\n",
			 s, PQ.ave_Z_V[s], PQ.ave_Z_V_err[s], max_nonzero(PQ.Z_V[s], N_K));
		}
		printf("  gxy^2l  : %6.2lf +- %.2lf  (max %d)\n",
		 PQ.ave_Z_g, PQ.ave_Z_g_err, max_nonzero(PQ.Z_g, N_K));
	}
	
	average_sub(PQ.free_energy, PQ.free_energy_err, n_bin, PQ.ave_sign);
	
	// constant term due to Canonical transformation
	// prm.ef[s] += real(D0_z_omega[0]) / 8.;
	PQ.free_energy += real(D0_omega[1][0]) / 8.;  // sp
	
// 	char str[32];
// 	sprintf(str, "\n V^2k up :");  // 10 letters
// 	analys_stat(PQ.Z_V[0], PQ.ave_Z_V[0], N_K, str);
// 	sprintf(str, " V^2k dw :");
// 	analys_stat(PQ.Z_V[1], PQ.ave_Z_V[1], N_K, str);
// 	sprintf(str, " gxy^2l  :");
// 	analys_stat(PQ.Z_g, PQ.ave_Z_g, N_K, str);
}

static void average_sp(int n_bin)
{
	average_sub(PQ.occup_tot, PQ.occup_tot_err, n_bin, PQ.ave_sign);
	average_sub(PQ.occup_dif, PQ.occup_dif_err, n_bin, PQ.ave_sign);
	average_sub(PQ.occup_dbl, PQ.occup_dbl_err, n_bin, PQ.ave_sign);
	
	for(int s=0; s<2; s++){
		average_sub(PQ.occup_ud[s], PQ.occup_ud_err[s], n_bin, PQ.ave_sign);
		
		for(int i=1; i<N_TAU; i++){
			average_sub(SP[s].Gf_tau[i], SP[s].Gf_tau_err[i], n_bin, PQ.ave_sign);
		}
		
		SP[s].Gf_tau[N_TAU] = -PQ.occup_ud[s];
		SP[s].Gf_tau_err[N_TAU] = PQ.occup_ud_err[s];
		
		if( UINF ){
			SP[s].Gf_tau[0] = PQ.occup_tot - 1.0;
			SP[s].Gf_tau_err[0] = PQ.occup_tot_err;
		}
		else{
			SP[s].Gf_tau[0] = -SP[s].Gf_tau[N_TAU] -1.0;
			SP[s].Gf_tau_err[0] = SP[s].Gf_tau_err[N_TAU];
		}
		
		SP[s].jump = - SP[s].Gf_tau[0] - SP[s].Gf_tau[N_TAU];
	}
	
	average_sub(PQ.m_z_sqr, PQ.m_z_sqr_err, n_bin, PQ.ave_sign);
	average_sub(PQ.m_z_pw4, PQ.m_z_pw4_err, n_bin, PQ.ave_sign);
	
	// FFT tau -> omega
	for(int s=0; s<2; s++){
// 		SP[s].jump = 1.0 - PQ.occup_ud[(s+1)%2];
		fft_fermion_radix2_tau2omega(SP[s].Gf_tau, SP[s].Gf_omega, prm.beta, N_TAU, SP[s].jump);
	}
	
	// self-energy (by Dyson equation)
	for(int s=0; s<2; s++){
		for(int i=0; i<N_TAU/2; i++){
			complex<double> i_omega_f = IMAG * (double)(2*i+1) * M_PI / prm.beta;
			
			SP[s].self_f_1[i] = i_omega_f - prm_ef_unrenorm[s] - Delta_omega[s][i] - 1.0 / SP[s].Gf_omega[i];
			SP[s].self_f_2[i] = i_omega_f - ( prm_ef_unrenorm[s] + Delta_omega[s][i] + 1.0 / SP[s].Gf_omega[i] ) * SP[s].jump;
			
// 			#if PHONON
// 			SP[s].self_f[i] += prm.g * prm.g / prm.w_ph;
// 			#endif // PHONON
		}
	}
	
	// self-energy (by improved estimators)
	for(int s=0; s<2; s++){
		for(int i=1; i<N_TAU; i++){
			average_sub(SP[s].G_self[i], n_bin, PQ.ave_sign);
		}
// 		SP[s].G_self[0] = (PQ.occup_dbl - PQ.occup_ud[(s+1)%2]) * prm_U_unrenorm;
// 		SP[s].G_self[N_TAU] = -PQ.occup_dbl * prm_U_unrenorm;
		SP[s].G_self[0] = 2.* SP[s].G_self[1] - SP[s].G_self[2];
		SP[s].G_self[N_TAU] = 2.* SP[s].G_self[N_TAU-1] - SP[s].G_self[N_TAU-2];
		
		double jump =  - SP[s].G_self[0] - SP[s].G_self[N_TAU];
		fft_fermion_radix2_tau2omega(SP[s].G_self, SP[s].self_f, prm.beta, N_TAU, jump);
		
		for(int i=0; i<N_TAU/2; i++){
// 			SP[s].self_f[i] *= prm.U;
			SP[s].self_f[i] /= SP[s].Gf_omega[i];
		}
	}
	
}

static void average_pm(int n_bin)
{
	//
	// chi_{+-}
	//
	
	// static
	average_sub(PQ.stat_suscep_pm, PQ.stat_suscep_pm_err, n_bin, PQ.ave_sign);
	
	// dynamical
// 	for(int i=1; i<N_TAU; i++){
	for(int i=0; i<=N_TAU; i++){
		average_sub(TP->chi_pm_tau[i], TP->chi_pm_tau_err[i], n_bin, PQ.ave_sign);
	}
// 	TP->chi_pm_tau[0] = PQ.occup_ud[0] / 2.;
// 	TP->chi_pm_tau_err[0] = PQ.occup_ud_err[0] / 2.;
// 	TP->chi_pm_tau[N_TAU] = PQ.occup_ud[1] / 2.;
// 	TP->chi_pm_tau_err[N_TAU] = PQ.occup_ud_err[1] / 2.;
	
	// FFT tau -> omega
	double jump = TP->chi_pm_tau[N_TAU] - TP->chi_pm_tau[0];
	fft_boson_radix2_tau2omega(TP->chi_pm_tau, TP->chi_pm_omega, prm.beta, N_TAU, jump);
	
	average_sub(PQ.m_xy_sqr, PQ.m_xy_sqr_err, n_bin, PQ.ave_sign);
	average_sub(PQ.m_xy_pw4, PQ.m_xy_pw4_err, n_bin, PQ.ave_sign);
	
	average_sub(PQ.m_rt_sqr, PQ.m_rt_sqr_err, n_bin, PQ.ave_sign);
	average_sub(PQ.m_rt_pw4, PQ.m_rt_pw4_err, n_bin, PQ.ave_sign);
	
	// symmetrize
// 	double chi_pm_sym[N_TAU+1];
// 	for(int i=0; i<=N_TAU; i++){
// 		chi_pm_sym[i] = chi_pm_sym[N_TAU-i] = (TP->chi_pm_tau[i] + TP->chi_pm_tau[N_TAU-i]) / 2.;
// 	}
// 	fft_boson_radix2_tau2omega(chi_pm_sym, TP->chi_pm_omega, prm.beta, N_TAU);
}

static void average_tp(int n_bin)
{
	// static
	average_sub(PQ.stat_suscep_sp, PQ.stat_suscep_sp_err, n_bin, PQ.ave_sign);
	average_sub(PQ.stat_suscep_ch, PQ.stat_suscep_ch_err, n_bin, PQ.ave_sign);
	
// 	PQ.stat_suscep_sp -= prm.beta * PQ.occup_dif * PQ.occup_dif / 4.;  // < S_z^2 >
// 	PQ.stat_suscep_ch -= prm.beta * PQ.occup_tot * PQ.occup_tot;
	
	// dynamical
	for(int i=0; i<=N_TP; i++){
		average_sub(TP->chi_sp_tau[i], TP->chi_sp_tau_err[i], n_bin, PQ.ave_sign);
		average_sub(TP->chi_ch_tau[i], TP->chi_ch_tau_err[i], n_bin, PQ.ave_sign);
		
// 		TP->chi_sp_tau[i] -= PQ.occup_dif * PQ.occup_dif / 4.;  // < S_z^2 >
// 		TP->chi_ch_tau[i] -= PQ.occup_tot * PQ.occup_tot;
	}
	
	// FFT tau -> omega
	fft_chi_after_interp(TP->chi_sp_tau, TP->chi_sp_omega);
	fft_chi_after_interp(TP->chi_ch_tau, TP->chi_ch_omega);
	
	// effective moment  (analytical continuation from n>0 to n=0)
	int n_pade = 100;
	if( n_pade > N_TP2-1 )  n_pade = N_TP2-1;
	complex<double> chi_sp_pade[n_pade], i_omega_b[n_pade];
	
	for(int i=0; i<n_pade; i++)  i_omega_b[i] = IMAG * (double)(2*(i+1)) * M_PI / prm.beta;
	
	pade_init(i_omega_b, TP->chi_sp_omega+1, chi_sp_pade, n_pade);
	complex<double> w = 0;
	PQ.stat_suscep_reg = real(pade_calc(i_omega_b, chi_sp_pade, w, n_pade));
	PQ.eff_mom = (PQ.stat_suscep_sp - PQ.stat_suscep_reg) / prm.beta;
	
}

static void average_vx(int n_bin)
{
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						average_sub(VX->Gfour_lo[s1][s2][iw1][iw2][nu], n_bin, PQ.ave_sign);
						average_sub(VX->gamma_lo[s1][s2][iw1][iw2][nu], n_bin, PQ.ave_sign);
					}
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					average_sub(VX->Gfour_tr[s][iw1][iw2][nu], n_bin, PQ.ave_sign);
					average_sub(VX->gamma_tr[s][iw1][iw2][nu], n_bin, PQ.ave_sign);
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			average_sub(VX->G[s][iw], n_bin, PQ.ave_sign);
			average_sub(VX->self[s][iw], n_bin, PQ.ave_sign);
			VX->self[s][iw] /= VX->G[s][iw];
		}
	}
	
	// spin symmetric
// 	complex<double> gf[2*N_VX1+N_VX2] = {0};
// 	complex<double> self[2*N_VX1+N_VX2] = {0};
// 	for(int s=0; s<2; s++){
// 		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
// 			gf[iw] += VX->G[s][iw] / 2.;
// 			self[iw] += VX->self[s][iw] / 2.;
// 		}
// 	}
	
	complex<double> chi0_lo[2][2*N_VX1][N_VX2];
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1; iw++){
			for(int nu=0; nu<N_VX2; nu++){
				chi0_lo[s][iw][nu] = - VX->G[s][iw] * VX->G[s][iw+nu];
			}
		}
	}
	complex<double> chi0_tr[2][2*N_VX1][N_VX2];
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1; iw++){
			for(int nu=0; nu<N_VX2; nu++){
				chi0_tr[s][iw][nu] = - VX->G[s][iw] * VX->G[(s+1)%2][iw+nu];
			}
		}
	}
	
	// vertex lo (from chi)
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
						VX->gamma2_lo[s1][s2][iw1][iw2][nu] = VX->Gfour_lo[s1][s2][iw1][iw2][nu];
					}
					// nu=0
					VX->gamma2_lo[s1][s2][iw1][iw2][0] -= VX->G[s1][iw1] * VX->G[s2][iw2];
				}
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int nu=0; nu<N_VX2; nu++){
				// iw1 == iw2,  s1 == s2
				VX->gamma2_lo[s1][s1][iw1][iw1][nu] -= chi0_lo[s1][iw1][nu];
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
// 						VX->gamma2_lo[s1][s2][iw1][iw2][nu]
// 						 /= -(chi0_lo[s1][iw1][nu] * chi0_lo[s2][iw2][nu]);
						VX->gamma2_lo[s1][s2][iw1][iw2][nu]
						 *= -prm.beta / (chi0_lo[s1][iw1][nu] * chi0_lo[s2][iw2][nu]);
					}
				}
			}
		}
	}
	
	// vertex tr (from chi)
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					VX->gamma2_tr[s][iw1][iw2][nu] = VX->Gfour_tr[s][iw1][iw2][nu];
				}
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int nu=0; nu<N_VX2; nu++){
				// iw1 == iw2
				VX->gamma2_tr[s][iw1][iw1][nu] -= chi0_tr[s][iw1][nu];
			}
		}
	}
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
// 					VX->gamma2_tr[s][iw1][iw2][nu]
// 					 /= -(chi0_tr[s][iw1][nu] * chi0_tr[s][iw2][nu]);
					VX->gamma2_tr[s][iw1][iw2][nu]
					 *= -prm.beta / (chi0_tr[s][iw1][nu] * chi0_tr[s][iw2][nu]);
				}
			}
		}
	}
	
	// vertex lo (direct evaluation)
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw1=0; iw1<2*N_VX1; iw1++){
				for(int iw2=0; iw2<2*N_VX1; iw2++){
					for(int nu=0; nu<N_VX2; nu++){
// 						VX->gamma_lo[s1][s2][iw1][iw2][nu] -= VX->self[s1][iw1] * VX->Gfour_lo[s1][s2][iw1][iw2][nu];
// 						VX->gamma_lo[s1][s2][iw1][iw2][nu] *= prm.beta / (VX->G[s1][iw1+nu] * chi0_lo[s2][iw2][nu]);
						VX->gamma_lo[s1][s2][iw1][iw2][nu] -= VX->self[s2][iw2+nu] * VX->Gfour_lo[s1][s2][iw1][iw2][nu];
						VX->gamma_lo[s1][s2][iw1][iw2][nu] *= prm.beta / (VX->G[s2][iw2] * chi0_lo[s1][iw1][nu]);
					}
				}
			}
		}
	}
	
	// vertex tr (direct evaluation)
	for(int s=0; s<2; s++){
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
// 					VX->gamma_tr[s][iw1][iw2][nu] -= VX->self[s][iw1] * VX->Gfour_tr[s][iw1][iw2][nu];
// 					VX->gamma_tr[s][iw1][iw2][nu] *= prm.beta / (VX->G[(s+1)%2][iw1+nu] * chi0_tr[s][iw2][nu]);
					VX->gamma_tr[s][iw1][iw2][nu] -= VX->self[(s+1)%2][iw2+nu] * VX->Gfour_tr[s][iw1][iw2][nu];
					VX->gamma_tr[s][iw1][iw2][nu] *= prm.beta / (VX->G[s][iw2] * chi0_tr[s][iw1][nu]);
				}
			}
		}
	}
	
	// susceptibility
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int nu=0; nu<N_VX2; nu++){
				average_sub(VX->suscep_lo[s1][s2][nu], n_bin, PQ.ave_sign);
			}
			VX->suscep_lo[s1][s2][0] -= PQ.occup_ud[s1] * PQ.occup_ud[s2] * prm.beta;
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					average_sub(VX->Gthree_lo[s1][s2][iw][nu], n_bin, PQ.ave_sign);
				}
				VX->Gthree_lo[s1][s2][iw][0] -= VX->G[s1][iw] * PQ.occup_ud[s2] * prm.beta;
			}
		}
	}
	
	// 3-point vertex
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					VX->lambda_lo[s1][s2][iw][nu] = VX->Gthree_lo[s1][s2][iw][nu];
// 					VX->lambda_lo[s1][s2][iw][nu] /= -VX->G[s1][iw] * VX->G[s1][iw+nu];
					VX->lambda_lo[s1][s2][iw][nu] /= chi0_lo[s1][iw][nu];
					
// 					if( s1 == s2 ){
// 						VX->lambda_lo[s1][s1][iw][nu] -= 1.0;
// 					}
				}
			}
		}
	}
	for(int s1=0; s1<2; s1++){
		for(int s2=0; s2<2; s2++){
			for(int iw=0; iw<2*N_VX1; iw++){
				for(int nu=0; nu<N_VX2; nu++){
					VX->lambda_lo_test[s1][s2][iw][nu] = 0;
					for(int iw2=0; iw2<2*N_VX1; iw2++){
// 						VX->lambda_lo_test[s1][s2][iw][nu] += VX->gamma_lo[s1][s2][iw][iw2][nu] * VX->G[s2][iw2] * VX->G[s2][iw2+nu];
						VX->lambda_lo_test[s1][s2][iw][nu] += VX->gamma_lo[s1][s2][iw][iw2][nu] * (-chi0_lo[s2][iw2][nu]);
					}
					VX->lambda_lo_test[s1][s2][iw][nu] /= prm.beta;
				}
			}
		}
	}
	for(int nu=0; nu<N_VX2; nu++){
		complex<double> denom = 1./ (VX->suscep_lo[0][0][nu] * VX->suscep_lo[1][1][nu] - VX->suscep_lo[0][1][nu] * VX->suscep_lo[1][0][nu]);
		complex<double> suscep_inv[2][2] = {
			{ VX->suscep_lo[1][1][nu] * denom, -VX->suscep_lo[0][1][nu] * denom},
			{-VX->suscep_lo[1][0][nu] * denom,  VX->suscep_lo[0][0][nu] * denom},
		};
		
		for(int iw=0; iw<2*N_VX1; iw++){
			complex<double> temp[2][2] = {
				{VX->lambda_lo[0][0][iw][nu], VX->lambda_lo[0][1][iw][nu]},
				{VX->lambda_lo[1][0][iw][nu], VX->lambda_lo[1][1][iw][nu]},
			};
			for(int s1=0; s1<2; s1++){
				for(int s2=0; s2<2; s2++){
					VX->lambda_lo[s1][s2][iw][nu] = 0;
					for(int s3=0; s3<2; s3++){
						VX->lambda_lo[s1][s2][iw][nu] += temp[s1][s3] * suscep_inv[s3][s2];
					}
				}
			}
		}
		for(int iw=0; iw<2*N_VX1; iw++){
			complex<double> temp[2][2] = {
				{VX->lambda_lo_test[0][0][iw][nu], VX->lambda_lo_test[0][1][iw][nu]},
				{VX->lambda_lo_test[1][0][iw][nu], VX->lambda_lo_test[1][1][iw][nu]},
			};
			for(int s1=0; s1<2; s1++){
				for(int s2=0; s2<2; s2++){
					VX->lambda_lo_test[s1][s2][iw][nu] = 0;
					for(int s3=0; s3<2; s3++){
						VX->lambda_lo_test[s1][s2][iw][nu] += temp[s1][s3] * suscep_inv[s3][s2];
					}
				}
			}
		}
	}
	
}
/*
static void average_vx(int n_bin)
{
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				average_sub(VX->chi_ch[iw1][iw2][nu], n_bin, PQ.ave_sign);
				average_sub(VX->chi_sp[iw1][iw2][nu], n_bin, PQ.ave_sign);
				average_sub(VX->vx_ch[iw1][iw2][nu], n_bin, PQ.ave_sign);
				average_sub(VX->vx_sp[iw1][iw2][nu], n_bin, PQ.ave_sign);
			}
		}
	}
	
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			average_sub(VX->G[s][iw], n_bin, PQ.ave_sign);
			average_sub(VX->self[s][iw], n_bin, PQ.ave_sign);
			VX->self[s][iw] /= VX->G[s][iw];
		}
	}
	
	// spin symmetric
	complex<double> gf[2*N_VX1+N_VX2] = {0};
	complex<double> self[2*N_VX1+N_VX2] = {0};
	for(int s=0; s<2; s++){
		for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
			gf[iw] += VX->G[s][iw] / 2.;
			self[iw] += VX->self[s][iw] / 2.;
		}
	}
	
	complex<double> chi0[2*N_VX1][N_VX2];
	for(int iw=0; iw<2*N_VX1; iw++){
		for(int nu=0; nu<N_VX2; nu++){
			chi0[iw][nu] = - gf[iw] * gf[iw+nu];
		}
	}
	
	// vertex
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				VX->vx2_ch[iw1][iw2][nu] = VX->chi_ch[iw1][iw2][nu];
				VX->vx2_sp[iw1][iw2][nu] = VX->chi_sp[iw1][iw2][nu];
			}
			
			// nu=0
// 			VX->vx2_ch[iw1][iw2][0] -= 2.* gf[iw1] * gf[iw2];
			VX->vx2_ch[iw1][iw2][0] -= 2.* gf[iw1] * gf[iw2] * prm.beta;
		}
		
		for(int nu=0; nu<N_VX2; nu++){
			// iw1 == iw2
// 			VX->vx2_ch[iw1][iw1][nu] -= chi0[iw1][nu];
// 			VX->vx2_sp[iw1][iw1][nu] -= chi0[iw1][nu];
			VX->vx2_ch[iw1][iw1][nu] -= chi0[iw1][nu] * prm.beta;
			VX->vx2_sp[iw1][iw1][nu] -= chi0[iw1][nu] * prm.beta;
		}
	}
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				complex<double> temp = chi0[iw1][nu] * chi0[iw2][nu];
				VX->vx2_ch[iw1][iw2][nu] /= -temp;  // *****sign
				VX->vx2_sp[iw1][iw2][nu] /= -temp;  // *****sign
			}
		}
	}
	
	// vertex
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				VX->vx_ch[iw1][iw2][nu] -= self[iw1+nu] * VX->chi_ch[iw1][iw2][nu];
				VX->vx_sp[iw1][iw2][nu] -= self[iw1+nu] * VX->chi_sp[iw1][iw2][nu];
				VX->vx_ch[iw1][iw2][nu] /= gf[iw1] * chi0[iw2][nu];
				VX->vx_sp[iw1][iw2][nu] /= gf[iw1] * chi0[iw2][nu];
			}
		}
	}
}
*/

static void func_average0(int n_bin)
{
	average_stat(n_bin);
}
static void func_average1(int n_bin)
{
	average_stat(n_bin);
	average_sp(n_bin);
	average_pm(n_bin);
}
static void func_average2(int n_bin)
{
	average_stat(n_bin);
	average_sp(n_bin);
	average_pm(n_bin);
	average_tp(n_bin);
}
static void func_average3(int n_bin)
{
	average_stat(n_bin);
	average_sp(n_bin);
	average_pm(n_bin);
	average_tp(n_bin);
	average_vx(n_bin);
}

//============================================================================
//
// MPI REDUCE
//

// return time
static double mpi_reduce_bin(int i_measure)
{
	double start=0, end=0;  // measure time
	
	#if QMC_MPI
	
// 	static struct phys_quant_bin{
// 		long int ave_sign;
// 		unsigned long n_k[2][N_K];  // V
// 		unsigned long n_l[N_K];  // g
// 		double Gf[2][N_TAU+1];
// 		double self[2][N_TAU+1];
// 		double occup[4];
// 		double m_z_sqr, m_z_pw4;
// 		double m_xy_sqr, m_xy_pw4;
// 		double m_rt_sqr, m_rt_pw4;
// 		double chi_pm[N_TAU+1];
// 		double chi[4][4][N_TP+1];
// 		complex<double> vx_Gf[2][2*N_VX1+N_VX2];
// 		complex<double> vx_self[2][2*N_VX1+N_VX2];
// 		complex<double> vx_G4_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];
// 		complex<double> vx_gamma_lo[2][2][2*N_VX1][2*N_VX1][N_VX2];
// 		complex<double> vx_G4_tr[2][2*N_VX1][2*N_VX1][N_VX2];
// 		complex<double> vx_gamma_tr[2][2*N_VX1][2*N_VX1][N_VX2];
// 		complex<double> vx_suscep[2][2][N_VX2];
// 		complex<double> vx_G3[2][2][2*N_VX1][N_VX2];
// 	} *B, *B_TOT;
	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	
	MPI_Reduce(&B->ave_sign, &B_TOT->ave_sign, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&B->n_k[0][0], &B_TOT->n_k[0][0], 2*N_K, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&B->n_l[0], &B_TOT->n_l[0], N_K, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(i_measure>0){
		MPI_Reduce(&B->Gf[0][0], &B_TOT->Gf[0][0], 2*(N_TAU+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->self[0][0], &B_TOT->self[0][0], 2*(N_TAU+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->occup[0], &B_TOT->occup[0], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->m_z_sqr, &B_TOT->m_z_sqr, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->chi_pm[0], &B_TOT->chi_pm[0], N_TAU+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	
	if( i_measure>1 ){
		MPI_Reduce(&B->chi[0][0][0], &B_TOT->chi[0][0][0], 4*4*(N_TP+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
// 		#if CHI_TR
// 		MPI_Reduce(&B->chi_tr1[0], &B_TOT->chi_tr1[0], (2*N_TP2+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		MPI_Reduce(&B->chi_tr2[0], &B_TOT->chi_tr2[0], (2*N_TP2+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		
// 		MPI_Reduce(&B->chi_omega[0], &B_TOT->chi_omega[0], (2*N_TP_W), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		#endif // CHI_TR
	}
	if( i_measure>2 ){
		MPI_Reduce(&B->vx_Gf[0][0], &B_TOT->vx_Gf[0][0], 4*(2*N_VX1+N_VX2), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->vx_self[0][0], &B_TOT->vx_self[0][0], 4*(2*N_VX1+N_VX2), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		// 2*2* (2*N_VX1)^2 * N_VX2 * 2 = 32*N_VX1*N_VX1*N_VX2
		MPI_Reduce(&B->vx_G4_lo[0][0][0][0][0], &B_TOT->vx_G4_lo[0][0][0][0][0], 32*N_VX1*N_VX1*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->vx_gamma_lo[0][0][0][0][0], &B_TOT->vx_gamma_lo[0][0][0][0][0], 32*N_VX1*N_VX1*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		// 2* (2*N_VX1)^2 * N_VX2 * 2 = 16*N_VX1*N_VX1*N_VX2
		MPI_Reduce(&B->vx_G4_tr[0][0][0][0], &B_TOT->vx_G4_tr[0][0][0][0], 16*N_VX1*N_VX1*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->vx_gamma_tr[0][0][0][0], &B_TOT->vx_gamma_tr[0][0][0][0], 16*N_VX1*N_VX1*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		// 2*2 * N_VX2 * 2 = 8 * N_VX2
		MPI_Reduce(&B->vx_suscep[0][0][0], &B_TOT->vx_suscep[0][0][0], 8*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		// 2*2* (2*N_VX1) * N_VX2 * 2 = 16*N_VX1*N_VX2
		MPI_Reduce(&B->vx_G3[0][0][0][0], &B_TOT->vx_G3[0][0][0][0], 16*N_VX1*N_VX2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	
	if(my_rank==0)  *B = *B_TOT;
	
	end = MPI_Wtime();
	
// 	if(my_rank==0)  print_time_mpi(end - start);
	
	#endif // QMC_MPI
	
	return(end - start);
	
}

static double mpi_reduce_accept()
{
	double start=0, end=0;  // measure time
	
	#if QMC_MPI
	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	
	//static struct mc_accept_sub{
	//	unsigned long n_accept_seg, n_reject_seg;
	//	unsigned long n_accept_spin, n_reject_spin;
	//	unsigned long n_accept_boson, n_reject_boson;
	//	unsigned long n_accept_shift, n_reject_shift;
	//	unsigned long n_accept_segsp, n_reject_segsp;
	//} ACCPT_sub;
	
	struct mc_accept_sub ACCPT_sub_tot;
	MPI_Reduce(&ACCPT_sub.n_accept_seg, &ACCPT_sub_tot.n_accept_seg, 10, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(my_rank==0){
		ACCPT_sub = ACCPT_sub_tot;
	}
	
	end = MPI_Wtime();
	
// 	if(my_rank==0)  print_time_mpi(end - start);
	
	#endif // QMC_MPI
	
	return(end - start);
}

//============================================================================
//
// UPDATES OF MC CONFIGURATIONS
//
//============================================================================


static void array_insert(double *arr, int N, double val, int pos)
{
	for(int i=N; i>pos; i--)  arr[i] = arr[i-1];
	arr[pos] = val;
}
static void array_insert(int *arr, int N, int val, int pos)
{
	for(int i=N; i>pos; i--)  arr[i] = arr[i-1];
	arr[pos] = val;
}
// if pos1 == pos2,  the order pos1 < pos2 is used
static void array_insert2(double *arr, int N, double val1, int pos1, double val2, int pos2)
{
	if( pos1 > pos2 ){
		array_insert(arr, N, val1, pos1);
		array_insert(arr, N+1, val2, pos2);
	}
	else{
		array_insert(arr, N, val2, pos2);
		array_insert(arr, N+1, val1, pos1);
	}
}
static void array_insert2(int *arr, int N, int val1, int pos1, int val2, int pos2)
{
	if( pos1 > pos2 ){
		array_insert(arr, N, val1, pos1);
		array_insert(arr, N+1, val2, pos2);
	}
	else{
		array_insert(arr, N, val2, pos2);
		array_insert(arr, N+1, val1, pos1);
	}
}

static void array_elim(double *arr, int N, int pos)
{
	for(int i=pos; i<N-1; i++)  arr[i] = arr[i+1];
}
static void array_elim(int *arr, int N, int pos)
{
	for(int i=pos; i<N-1; i++)  arr[i] = arr[i+1];
}
static void array_elim2(double *arr, int N, int pos1, int pos2)
{
	if( pos1 > pos2 ){
		array_elim(arr, N, pos1);
		array_elim(arr, N-1, pos2);
	}
	else if( pos1 < pos2 ){
		array_elim(arr, N, pos2);
		array_elim(arr, N-1, pos1);
	}
	else{
		printf("*** error array_elim2\n");
		exit(0);
	}
}
static void array_elim2(int *arr, int N, int pos1, int pos2)
{
	if( pos1 > pos2 ){
		array_elim(arr, N, pos1);
		array_elim(arr, N-1, pos2);
	}
	else if( pos1 < pos2 ){
		array_elim(arr, N, pos2);
		array_elim(arr, N-1, pos1);
	}
	else{
		printf("*** error array_elim2\n");
		exit(0);
	}
}

// static void array_rotate_upward(double *arr, int N, int n)
// {
// 	double temp[n];
// 	for(int i=0; i<n; i++)  temp[i] = arr[N-n+i];
// 	for(int i=N-1-n; i>=0; i--)  arr[i] = arr[i+n];
// 	for(int i=0; i<n; i++)  arr[i] = temp[i];
// }

static double length_occupied(int sigma, double tau1, double tau2)
{
	double len = 0;
	
	if( tau1 > tau2 ){
		len += length_occupied(sigma, 0.0, tau2);
		len += length_occupied(sigma, tau1, prm.beta);
	}
	else if( X.k ){
		int i = tau_order(X.tau, X.k, tau1);
		
		if( X.tau[i] > tau2 ){
			if( occupied(X.state[i], sigma) )  len += tau2 - tau1;
		}
		else{
			if( occupied(X.state[i], sigma) )  len += X.tau[i] - tau1;
			
			while( X.tau[++i] < tau2 ){
				if( occupied(X.state[i], sigma) )  len += X.tau[i] - X.tau[i-1];
			}
// 			while( X.tau[i+1] < tau2 ){
// 				if( occupied(X.state[i+1], sigma) )  len += X.tau[i+1] - X.tau[i];
// 				i++;
// 			}
			
			if( occupied(X.state[i], sigma) )  len += tau2 - X.tau[i-1];
		}
	}
	else{
		if( occupied(X.state[0], sigma) )  len += tau2 - tau1;
	}
	
	return len;
}
// static double length_occupied(int sigma, double tau1, double tau2)
// {
// 	double len = 0;
// 	
// 	if( tau1 > tau2 ){
// 		len += length_occupied(sigma, 0.0, tau2);
// 		len += length_occupied(sigma, tau1, prm.beta);
// 	}
// 	else{
// 		if( X.k ){
// 			int i_tau1 = tau_order(X.tau, X.k, tau1);
// 			int i_tau2 = tau_order(X.tau, X.k, tau2);
// 			
// 			if( occupied(X.state[i_tau1], sigma) )  len += X.tau[i_tau1] - tau1;
// 			if( occupied(X.state[i_tau2], sigma) )  len += tau2 - X.tau[i_tau2];
// 			len += length_occupied(sigma, i_tau1, i_tau2-1);
// 		}
// 		else{
// 			if( occupied(X.state[0], sigma) )  len += tau2 - tau1;
// 		}
// 	}
// 	return len;
// }
// static double length_occupied(int sigma, int i_tau1, int i_tau2)
// {
// 	double len = 0;
// 	if( i_tau1 > i_tau2 ){
// 		len += length_occupied(sigma, 0.0, i_tau2);
// 		len += length_occupied(sigma, i_tau1, prm.beta);
// 	}
// 	else{
// 		for(int i=i_tau1; i<i_tau2; i++){
// 			if( occupied(X.state[i], sigma) )  len += X.tau[i+1] - X.tau[i];
// 		}
// 	}
// 	return len;
// }

// return 1 if sigma segment can be created at position i
static int possible_create_seg(int i, int sigma, int anti)
{
	#if UINF
	if( occupied(X.state[i], (sigma+1)%2) )  return 0;
	#endif
	
	if( occupied(X.state[i], sigma) == anti )  return 1;
	return 0;
}

// static void swap_val(double &x1, double &x2)
// {
// 	double temp = x1;
// 	x1 = x2;
// 	x2 = temp;
// }
template <class T>
inline void swap_val(T& a, T& b)
{
	T tmp(a);
	a = b;
	b = tmp;
}

// anti = 0: segment     (f- f+)
// anti = 1: antisegment (f+ f-)
static void add_seg(int sigma, int anti)
{
	cond_op *F = &S[sigma];
	
	double tau_l, tau_r;
	int i_tau;
	
	//
	// choose tau_l & tau_r
	//
	tau_r = rand_tau_beta(prm.beta);
	i_tau = tau_order(X.tau, X.k, tau_r);
	
	#if TEST_MODE
	printf("\nadd_seg  (sigma=%d  anti=%d)\n", sigma, anti);
	printf("  i_tau=%d  tau_r=%.4lf\n", i_tau, tau_r);
	#endif
	
// 	int state = op2state[X.op[i_tau]];
// 	if( state != op2state[ F_OP[sigma][(anti+1)%2] ] )  return;
// 	// anti=0: f+ => empty
// 	// anti=1: f- => sigma
	
// 	if( FORBIDDEN(i_tau) )  return;
	if( ! possible_create_seg(i_tau, sigma, anti) )  return;
	
	double l_max = prm.beta;
	if( UINF ){
		if( X.k )  l_max = X.tau[i_tau] - tau_r;
	}
// 	if( X.k ){
	else if( F->k || Y.k ){
// 		int i = 0;
// 		for(; i<X.k; i++){
// 			if( ! possible_create_seg( (i_tau+i+1)%X.k, sigma, anti ) )  break;
// 		}
// 		l_max = X.tau[(i_tau+i)%X.k] - tau_r;
		
		int i = i_tau;
		while( possible_create_seg( (i+1)%X.k, sigma, anti ) ){  i++; }
		l_max = X.tau[i%X.k] - tau_r;
		
		if( l_max < 0 )  l_max += prm.beta;
	}
	double l_seg = rand_tau_l(l_max);
	tau_l = tau_r + l_seg;
	
	int flag_rotate = 0;
	if( tau_l >= prm.beta ){
		flag_rotate = 1;
		tau_l -= prm.beta;
	}
	int i_tau_l = tau_order(X.tau, X.k, tau_l);
	
	//
	// probability
	//
	double tau1;  // for f-annihilation operator
	double tau2;  // for f-creation operator
	
	if( !anti ){
		tau1 = tau_l;  tau2 = tau_r;  // f- f+
	} else {
		tau2 = tau_l;  tau1 = tau_r;  // f+ f-
	}
	
	double Delta_tau1[F->k], Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}
	
	double diag = G0_calc_interp(Delta[sigma], tau2, tau1);
	double lambda = add_lambda(*F, Delta_tau1, Delta_tau2, diag);
	
	double sign_b_fac = double(1 - 2*anti); // +1, -1
	double l_dbl = UINF ? 0 : length_occupied((sigma+1)%2, tau_r, tau_l);
	double b_fac = exp( -sign_b_fac * (l_seg * prm.ef[sigma] + l_dbl * prm.U) );
	
	double ln_z_fac = 0;
	for(int a=0; a<2; a++){
		if( flag_g[a] ){
			double s_fac1 = s_fac[a][F_OP[sigma][0]];  // tau1 f-
			double s_fac2 = s_fac[a][F_OP[sigma][1]];  // tau2 f+
			ln_z_fac += D0_calc_interp(K0[a], tau2, tau1) * s_fac1 * s_fac2;
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau1) * s_fac[a][X.op[i]] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau2) * s_fac[a][X.op[i]] * s_fac2;
			}
		}
	}
	double z_fac = exp( ln_z_fac );
	
	double prob = lambda * b_fac * z_fac * prm.beta * l_max / (double)(F->k+1);
	if( anti ) prob = -prob;
	if( flag_rotate )  prob = -prob;
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  l_max=%.5lf\n", l_max);
	printf("  lambda=%.5lf\n", lambda);
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	//
	// update
	//
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		
		add_mat_M(*F, F->k, F->k, Delta_tau1, Delta_tau2, lambda);
		F->tau1[F->k] = tau1;
		F->tau2[F->k] = tau2;
		F->k ++;
		
		int F_l = F_OP[sigma][anti];
		int F_r = F_OP[sigma][(anti+1)%2];
		array_insert(X.tau, X.k, tau_l, i_tau_l);
		array_insert(X.tau, X.k+1, tau_r, i_tau+flag_rotate);
// 		array_insert2(X.tau, X.k, tau_r, i_tau, tau_l, i_tau_l);
		array_insert(X.op, X.k, F_l, i_tau_l);
		array_insert(X.op, X.k+1, F_r, i_tau+flag_rotate);
// 		array_insert2(X.op, X.k, F_r, i_tau, F_l, i_tau_l);
		int state_l = X.state[i_tau_l]^(0x01 << sigma);
		int state_r = X.state[i_tau];
		
		int n_bw = i_tau_l - i_tau + flag_rotate * X.k;
		for(int i=0; i<n_bw; i++){
			X.state[(i_tau+i)%X.k] ^= (0x01 << sigma);
		}
		array_insert(X.state, X.k, state_l, i_tau_l);
		array_insert(X.state, X.k+1, state_r, i_tau+flag_rotate);
// 		array_insert2(X.state, X.k, state_r, i_tau, state_l, i_tau_l);
		X.k += 2;
		
// 		if( flag_rotate ){
// 		if( X.k == 2 && flag_rotate ){
// 			rotate_upward(X.tau, X.k);
// 			rotate_upward(X.op, X.k);
// 			rotate_upward(X.state, X.k);
// 		}
		
		X.tau[X.k] = X.tau[0] + prm.beta;
		X.tau[X.k+1] = X.tau[1] + prm.beta;
		X.op[X.k] = X.op[0];
		X.state[X.k] = X.state[0];
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}

static void rem_seg(int sigma, int anti)
{
	cond_op *F = &S[sigma];
	if( F->k == 0 )  return;
	
	int j_tau_l, j_tau_r;
	double tau_l, tau_r;
	double *F_tau_l, *F_tau_r;
	if( !anti){
		F_tau_l = F->tau1;  F_tau_r = F->tau2;  // f- f+
	} else {
		F_tau_l = F->tau2;  F_tau_r = F->tau1;  // f+ f-
	}
	
	//
	// choose tau_l & tau_r
	//
	j_tau_r = rand_int(F->k);
	tau_r = F_tau_r[j_tau_r];
	
	int i_tau = tau_position(X.tau, X.k, tau_r);
	
	#if TEST_MODE
	printf("\nrem_seg  (sigma=%d  anti=%d)\n", sigma, anti);
// 	printf(" j_tau2=%d\n", j_tau2);
	printf(" i_tau=%d\n", i_tau);
	#endif
	
// 	if( X.op[i_tau+1] != F_OP[sigma][anti] )  return;
	// anti=0: f-
	// anti=1: f+
	
	int i_tau_l = i_tau;
	while( X.op[(++i_tau_l)%X.k] != F_OP[sigma][anti] ){
		if( UINF )  return;
		if( X.op[i_tau_l%X.k] == S_PLUS || X.op[i_tau_l%X.k] == S_MINUS )  return;
	}
	
// 	tau_l = X.tau[i_tau+1];
// 	int flag_rotate = 0;
// 	if( i_tau+1 == X.k ){
// 		flag_rotate = 1;
// 		tau_l = X.tau[0];
// 	}
	int flag_rotate = 0;
	if( i_tau_l >= X.k ){
		flag_rotate = 1;
		i_tau_l -= X.k;
// 		tau_l = X.tau[0];
	}
	tau_l = X.tau[i_tau_l];
	
	j_tau_l = tau_position(F_tau_l, F->k, tau_l);
	
	//
	// probability
	//
// 	double l_max = X.tau[i_tau+2] - X.tau[i_tau];
	double l_max = prm.beta;
	if( UINF ){
		if( X.k > 2 )  l_max = X.tau[i_tau+2] - X.tau[i_tau];
	}
// 	if( X.k > 2 ){
	else if( F->k > 1 || Y.k ){
// 		int i = 0;
// 		for(; i<X.k; i++){
// 			if( ! possible_create_seg( (i_tau_l+i+2)%X.k, sigma, anti ) )  break;
// 		}
// 		l_max = X.tau[(i_tau_l+i+1)%X.k] - tau_r;
		
		int i = i_tau_l+1;
		while( possible_create_seg( (i+1)%X.k, sigma, anti ) ){  i++; }
		l_max = X.tau[i%X.k] - tau_r;
		
		if( l_max < 0 )  l_max += prm.beta;
	}
	
// 	double l_seg = X.tau[i_tau+1] - X.tau[i_tau];
	double l_seg = X.tau[i_tau_l] - X.tau[i_tau];
	if( l_seg < 0 )  l_seg += prm.beta;
	double l_dbl = UINF ? 0 : length_occupied((sigma+1)%2, X.tau[i_tau], X.tau[i_tau_l]);
	double sign_b_fac = double(1 - 2*anti); // +1, -1
	double b_fac = exp( sign_b_fac * (l_seg * prm.ef[sigma] + l_dbl * prm.U) );
	
	int j_tau1, j_tau2;
	if( !anti){
		j_tau1 = j_tau_l;  j_tau2 = j_tau_r;
	} else {
		j_tau2 = j_tau_l;  j_tau1 = j_tau_r;
	}
	double lambda = F->mat_M[j_tau1][j_tau2];
	
	double ln_z_fac = 0;
	for(int a=0; a<2; a++){
		if( flag_g[a] ){
			double s_fac1 = s_fac[a][F_OP[sigma][anti]];  // tau_l
			double s_fac2 = s_fac[a][F_OP[sigma][(anti+1)%2]];  // tau_r
			ln_z_fac -= D0_calc_interp(K0[a], tau_l, tau_r) * s_fac1 * s_fac2;  // subtract double counting
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_l) * s_fac[a][X.op[i]] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_r) * s_fac[a][X.op[i]] * s_fac2;
			}
		}
	}
	double z_fac = exp( -ln_z_fac );
	
	double prob = lambda * b_fac * z_fac * (double)(F->k) / (prm.beta * l_max);
	if( anti ) prob = -prob;
	if( flag_rotate )  prob = -prob;
// 	if( (j_tau1 + j_tau2) % 2 )  prob = -prob;
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  j_tau1=%d  j_tau2=%d\n", j_tau1, j_tau2);
	printf("  l_max=%.5lf\n", l_max);
	printf("  lambda=%.5lf\n", lambda);
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	//
	// update
	//
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		
		remove_mat_M(*F, j_tau1, j_tau2, lambda);
		array_elim(F->tau1, F->k, j_tau1);
		array_elim(F->tau2, F->k, j_tau2);
		F->k --;
		
		if( X.k > 2 ){
// 			int i_tau_2 = i_tau;
// 			if( flag_rotate )  i_tau_2 = 0;
// 			array_elim(X.tau, X.k, i_tau);
// 			array_elim(X.tau, X.k-1, i_tau_2);
// 			array_elim(X.op, X.k, i_tau);
// 			array_elim(X.op, X.k-1, i_tau_2);
// 			array_elim(X.state, X.k, i_tau);
// 			array_elim(X.state, X.k-1, i_tau_2);
			
			int n_bw = i_tau_l - i_tau - 1 + flag_rotate * X.k;
			for(int i=0; i<n_bw; i++){
				X.state[(i_tau+i+1)%X.k] ^= (0x01 << sigma);
			}
			
			array_elim2(X.tau, X.k, i_tau, i_tau_l);
			array_elim2(X.op, X.k, i_tau, i_tau_l);
			array_elim2(X.state, X.k, i_tau, i_tau_l);
			X.k -= 2;
			
			X.tau[X.k] = X.tau[0] + prm.beta;
			X.tau[X.k+1] = X.tau[1] + prm.beta;
			X.op[X.k] = X.op[0];
			X.state[X.k] = X.state[0];
		}
		else{
// // 			X.tau[0] = prm.beta;
// 			X.op[0] = F_OP[sigma][(anti+1)%2];
			// anti=0: f+ => empty
			// anti=1: f- => sigma
			X.state[0] = X.state[i_tau];
			X.k = 0;
		}
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}

static void add_spin(int anti)
{
	double tau_l, tau_r;
	int i_tau;
	
	//
	// choose tau_l & tau_r
	//
	tau_r = rand_tau_beta(prm.beta);
	i_tau = tau_order(X.tau, X.k, tau_r);
	
	#if TEST_MODE
	printf("\nadd_spin  (anti=%d)\n", anti);
	printf("  i_tau=%d  tau_r=%.4lf\n", i_tau, tau_r);
	#endif
	
// 	int state = op2state[X.op[i_tau]];
// 	if( state != op2state[ S_OP[(anti+1)%2] ] )  return;
	// anti=0: S+ => DW
	// anti=1: S- => UP
	int sigma = (anti+1)%2;
	if( X.state[i_tau] != 0x01<<sigma )  return;
	
	double l_max = X.tau[i_tau] - tau_r;
	if( X.k == 0 )  l_max = prm.beta;
	double l_seg = rand_tau_l(l_max);
	tau_l = tau_r + l_seg;
	
	int flag_rotate = 0;
	if( tau_l >= prm.beta ){
		flag_rotate = 1;
		tau_l -= prm.beta;
	}
	
	//
	// probability
	//
	double tau1;  // S-
	double tau2;  // S+
	if( !anti ){
		tau1 = tau_l;  tau2 = tau_r;  // S- S+
	} else {
		tau2 = tau_l;  tau1 = tau_r;  // S+ S-
	}
	
	double lambda = D0_calc_interp(D0_pm, tau1, tau2);
	
	double sign_b_fac = double(1 - 2*anti); // +1, -1
	double b_fac = exp( -sign_b_fac * l_seg * (prm.ef[0] - prm.ef[1]) );
	
	double ln_z_fac = 0;
	{
		int a = 1;  // only sp
		if( flag_g[a] ){
			double s_fac1 = s_fac[a][S_MINUS];  // tau1 S-
			double s_fac2 = s_fac[a][S_PLUS];  // tau2 S+
			ln_z_fac += D0_calc_interp(K0[a], tau2, tau1) * s_fac1 * s_fac2;
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau1) * s_fac[a][X.op[i]] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau2) * s_fac[a][X.op[i]] * s_fac2;
			}
		}
	}
	double z_fac = exp( ln_z_fac );
	
	double prob = -0.5 * lambda * b_fac * z_fac * prm.beta * l_max / (double)(Y.k+1);
	prob *= exp( WL.log_gk[Y.k] - WL.log_gk[Y.k+1] );
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  l_max=%.5lf\n", l_max);
	printf("  lambda=%.5lf\n", lambda);
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	//
	// update
	//
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_spin, ACCPT_sub.n_reject_spin) ){
		
		Y.tau1[Y.k] = tau1;
		Y.tau2[Y.k] = tau2;
		Y.k++;
		
		array_insert(X.tau, X.k, tau_l, i_tau);
		array_insert(X.tau, X.k+1, tau_r, i_tau);
		array_insert(X.op, X.k, S_OP[anti], i_tau);
		array_insert(X.op, X.k+1, S_OP[(anti+1)%2], i_tau);
		int state_l = X.state[i_tau]^0x03;  // UP <--> DW
		int state_r = X.state[i_tau];
		array_insert(X.state, X.k, state_l, i_tau);
		array_insert(X.state, X.k+1, state_r, i_tau);
		X.k += 2;
		
		if( flag_rotate ){
			rotate_upward(X.tau, X.k);
			rotate_upward(X.op, X.k);
			rotate_upward(X.state, X.k);
		}
		
		X.tau[X.k] = X.tau[0] + prm.beta;
		X.tau[X.k+1] = X.tau[1] + prm.beta;
		X.op[X.k] = X.op[0];
		X.state[X.k] = X.state[0];
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}

static void rem_spin(int anti)
{
	if( Y.k == 0 )  return;
	
	int j_tau_l, j_tau_r;
	double tau_l, tau_r;
	double *D_tau_l, *D_tau_r;
	if( !anti ){
		D_tau_l = Y.tau1;  D_tau_r = Y.tau2;  // S- S+
	} else {
		D_tau_l = Y.tau2;  D_tau_r = Y.tau1;  // S+ S-
	}
	
	//
	// choose tau_l & tau_r
	//
	j_tau_r = rand_int(Y.k);
	tau_r = D_tau_r[j_tau_r];
	
	int i_tau = tau_position(X.tau, X.k, tau_r);
	
	#if TEST_MODE
	printf("\nrem_spin  (anti=%d)\n", anti);
	printf(" j_tau_r=%d\n", j_tau_r);
	printf(" i_tau=%d\n", i_tau);
	#endif
	
	if( X.op[i_tau+1] != S_OP[anti] )  return;
	// anti=0: S-
	// anti=1: S+
// 	printf("ok1\n");
	
	tau_l = X.tau[i_tau+1];
	int flag_rotate = 0;
	if( i_tau+1 == X.k ){
		flag_rotate = 1;
		tau_l = X.tau[0];
	}
	
	j_tau_l = tau_position(D_tau_l, Y.k, tau_l);
	if( j_tau_l != j_tau_r )  return;
// 	printf("ok2\n");
	
	//
	// probability
	//
	double l_max = X.tau[i_tau+2] - X.tau[i_tau];
	double l_seg = X.tau[i_tau+1] - X.tau[i_tau];
	
	double tau1, tau2;
	if( !anti ){
		tau1 = tau_l;  tau2 = tau_r;  // S- S+
	} else {
		tau2 = tau_l;  tau1 = tau_r;  // S+ S-
	}
	double lambda = 1./ D0_calc_interp(D0_pm, tau1, tau2);  // symmetric
	
	double sign_b_fac = double(1 - 2*anti); // +1, -1
	double b_fac = exp( sign_b_fac * l_seg * (prm.ef[0] - prm.ef[1]) );
	
	double ln_z_fac = 0;
	{
		int a = 1;  // only sp
		if( flag_g[a] ){
			double s_fac1 = s_fac[a][S_OP[anti]];  // tau_l
			double s_fac2 = s_fac[a][S_OP[(anti+1)%2]];  // tau_r
			ln_z_fac -= D0_calc_interp(K0[a], tau_l, tau_r) * s_fac1 * s_fac2;  // subtract double counting
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_l) * s_fac[a][X.op[i]] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_r) * s_fac[a][X.op[i]] * s_fac2;
			}
		}
	}
	double z_fac = exp( -ln_z_fac );
	
	double prob = -2.* lambda * b_fac * z_fac * double(Y.k) / (prm.beta * l_max);
	prob *= exp( WL.log_gk[Y.k] - WL.log_gk[Y.k-1] );
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
// 	printf("  j_tau1=%d  j_tau2=%d\n", j_tau1, j_tau2);
	printf("  l_max=%.5lf\n", l_max);
	printf("  lambda=%.5lf\n", lambda);
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	//
	// update
	//
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_spin, ACCPT_sub.n_reject_spin) ){
		
		array_elim(D_tau_l, Y.k, j_tau_l);  // j_tau_l == j_tau_r
		array_elim(D_tau_r, Y.k, j_tau_r);
		Y.k --;
		
		if( X.k > 2 ){
			int i_tau_2 = i_tau;
			if( flag_rotate )  i_tau_2 = 0;
			array_elim(X.tau, X.k, i_tau);
			array_elim(X.tau, X.k-1, i_tau_2);
			array_elim(X.op, X.k, i_tau);
			array_elim(X.op, X.k-1, i_tau_2);
			array_elim(X.state, X.k, i_tau);
			array_elim(X.state, X.k-1, i_tau_2);
			X.k -= 2;
			
			X.tau[X.k] = X.tau[0] + prm.beta;
			X.tau[X.k+1] = X.tau[1] + prm.beta;
			X.op[X.k] = X.op[0];
			X.state[X.k] = X.state[0];
		}
		else{
// // 			X.tau[0] = prm.beta;
// 			X.op[0] = S_OP[(anti+1)%2];
			// anti=0: S+ => DW
			// anti=1: S- => UP
			X.state[0] = X.state[i_tau];
			X.k = 0;
		}
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}

static void spin2seg(int dbl_m, int dbl_p)
{
	if( Y.k == 0 )  return;
	
	int F_m_l = F_DW_CR, F_m_r = F_UP_AN;  // S- -> f+dw f-up
	int F_p_l = F_UP_CR, F_p_r = F_DW_AN;  // S+ -> f+up f-dw
	int state_m = EMP, state_p = EMP;
	if( dbl_m ){  swap_val(F_m_l, F_m_r);  state_m = DBL;}
	if( dbl_p ){  swap_val(F_p_l, F_p_r);  state_p = DBL;}
	
	//
	// choose tau_1 & tau_2
	//
	int j_tau = rand_int(Y.k);
// 	double tau_m = tau1[0] = Y.tau1[j_tau];  // S- -> f-up
// 	double tau_p = tau1[1] = Y.tau2[j_tau];  // S+ -> f-down
	double tau_m = Y.tau1[j_tau];  // S- -> f-up
	double tau_p = Y.tau2[j_tau];  // S+ -> f-down
	
	int i_tau_m = tau_position(X.tau, X.k, tau_m);
	int i_tau_p = tau_position(X.tau, X.k, tau_p);
	
	double l_max_m = X.tau[i_tau_m+1] - tau_m;
	double l_seg_m = rand_tau_l(l_max_m);
// 	tau2[1] = tau_m + l_seg_m;  // f+down
	double tau_m_ins = tau_m + l_seg_m;  // f+down
	
	double l_max_p = X.tau[i_tau_p+1] - tau_p;
	double l_seg_p = rand_tau_l(l_max_p);
// 	tau2[0] = tau_p + l_seg_p;  // f+up
	double tau_p_ins = tau_p + l_seg_p;  // f+up
	
	int flag_rotate = 0;
// 	if( tau2[1] >= prm.beta ){
	if( tau_m_ins >= prm.beta ){
		flag_rotate = 1;
		tau_m_ins -= prm.beta;
	}
// 	if( tau2[0] >= prm.beta ){
	if( tau_p_ins >= prm.beta ){
		flag_rotate = 1;
		tau_p_ins -= prm.beta;
	}
	
	double tau1[2];  // f-[sigma]
	double tau2[2];  // f+[sigma]
	tau2[1] = tau_m_ins;  tau1[0] = tau_m;  // S- -> f+dw f-up
	tau2[0] = tau_p_ins;  tau1[1] = tau_p;  // S+ -> f+up f-dw
	if( dbl_m )  swap_val(tau2[1], tau1[0]);
	if( dbl_p )  swap_val(tau2[0], tau1[1]);
	
	#if TEST_MODE
	printf("\nspin2seg  (dbl_m=%d  dbl_p=%d)\n", dbl_m, dbl_p);
	printf(" j_tau=%d\n", j_tau);
	printf(" tau_m  =%.4lf  tau_p  =%.4lf\n", tau_m, tau_p);
	printf(" tau1[0]=%.4lf  tau1[1]=%.4lf\n", tau1[0], tau1[1]);
	printf(" tau2[1]=%.4lf  tau2[0]=%.4lf\n", tau2[1], tau2[0]);
	#endif
	
	//
	// probability
	//
	double lambda_b = D0_calc_interp(D0_pm, tau_m, tau_p);
	
	// l_seg_m : DW -> EMP (if dbl_m==0),  DBL (if dbl_m==1)
	// l_seg_p : UP -> EMP (if dbl_m==0),  DBL (if dbl_m==1)
	double ef_m = ( dbl_m==0 ) ? -prm.ef[1] : prm.ef[0] + prm.U;
	double ef_p = ( dbl_p==0 ) ? -prm.ef[0] : prm.ef[1] + prm.U;
	double b_fac = exp( -l_seg_m * ef_m - l_seg_p * ef_p );
	
	
	int k_max = MAX( S[0].k, S[1].k );
	double Delta_tau1[2][k_max], Delta_tau2[2][k_max];
	double diag[2], lambda_f[2];
	for(int s=0; s<2; s++){
		for(int i=0; i<S[s].k; i++){
			Delta_tau1[s][i] = G0_calc_interp(Delta[s], S[s].tau2[i], tau1[s]);
			Delta_tau2[s][i] = G0_calc_interp(Delta[s], tau2[s], S[s].tau1[i]);
		}
		
		diag[s] = G0_calc_interp(Delta[s], tau2[s], tau1[s]);
		lambda_f[s] = add_lambda(S[s], Delta_tau1[s], Delta_tau2[s], diag[s]);
	}
	
// 	int F_m_l = F_DW_CR, F_m_r = F_UP_AN;  // S- -> f+dw f-up
// 	int F_p_l = F_UP_CR, F_p_r = F_DW_AN;  // S+ -> f+up f-dw
	double ln_z_fac = 0;
	for(int a=0; a<2; a++){  // replace (S-, S+) with (f-up, f-down)
		if( flag_g[a] ){
	// 		double s_fac1 = s_fac_sp[F_UP_AN] - s_fac_sp[S_MINUS];  // tau_m : S- -> f-up
	// 		double s_fac2 = s_fac_sp[F_DW_AN] - s_fac_sp[S_PLUS];   // tau_p : S+ -> f-down
			double s_fac1 = s_fac[a][F_m_r] - s_fac[a][S_MINUS];  // tau_m : S- -> f-up
			double s_fac2 = s_fac[a][F_p_r] - s_fac[a][S_PLUS];   // tau_p : S+ -> f-down
			ln_z_fac += D0_calc_interp(K0[a], tau_p, tau_m)
	// 		         * (s_fac_sp[F_UP_AN] * s_fac_sp[F_DW_AN] - s_fac_sp[S_MINUS] * s_fac_sp[S_PLUS]);
			         * (s_fac[a][F_m_r] * s_fac[a][F_p_r] - s_fac[a][S_MINUS] * s_fac[a][S_PLUS]);
			
			double s_fac_cav[X.k];
			for(int i=0; i<X.k; i++)  s_fac_cav[i] = s_fac[a][X.op[i]];
			s_fac_cav[i_tau_m] = s_fac_cav[i_tau_p] = 0;
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_m) * s_fac_cav[i] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_p) * s_fac_cav[i] * s_fac2;
			}
		}
	}
	for(int a=0; a<2; a++){  // add (f+up, f+down)
		if( flag_g[a] ){
	// 		double s_fac1 = s_fac_sp[F_UP_CR];  // tau2[0] f+up
	// 		double s_fac2 = s_fac_sp[F_DW_CR];  // tau2[1] f+down
			double s_fac1 = s_fac[a][F_p_l];  // tau2[0] f+up
			double s_fac2 = s_fac[a][F_m_l];  // tau2[1] f+down
			
			double s_fac_new[X.k];
			for(int i=0; i<X.k; i++)  s_fac_new[i] = s_fac[a][X.op[i]];
	// 		s_fac_new[i_tau_m] = s_fac_sp[F_UP_AN];
	// 		s_fac_new[i_tau_p] = s_fac_sp[F_DW_AN];
			s_fac_new[i_tau_m] = s_fac[a][F_m_r];
			s_fac_new[i_tau_p] = s_fac[a][F_p_r];
	// 		ln_z_fac += D0_calc_interp(K0_sp, tau2[0], tau2[1]) * s_fac1 * s_fac2;
			ln_z_fac += D0_calc_interp(K0[a], tau_p_ins, tau_m_ins) * s_fac1 * s_fac2;
			for(int i=0; i<X.k; i++){
	// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[0]) * s_fac_new[i] * s_fac1;
	// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[1]) * s_fac_new[i] * s_fac2;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_p_ins) * s_fac_new[i] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_m_ins) * s_fac_new[i] * s_fac2;
			}
		}
	}
	double z_fac = exp( ln_z_fac );
	
// 	double ln_z_fac = 0;
// 	{
// 		for(int i=0; i<X.k; i++){
// 			for(int j=i+1; j<X.k; j++){
// 				ln_z_fac -= D0_calc_interp(K0_sp, X.tau[i], X.tau[j]) * s_fac_sp[X.op[i]] * s_fac_sp[X.op[j]];
// 			}
// 		}
// 		
// 		double s_fac_new[X.k+2];
// 		double tau_new[X.k+2];
// 		for(int i=0; i<X.k; i++){
// 			s_fac_new[i] = s_fac_sp[X.op[i]];
// 			tau_new[i] = X.tau[i];
// 		}
// 		s_fac_new[i_tau_m] = s_fac_sp[F_UP_AN];
// 		s_fac_new[i_tau_p] = s_fac_sp[F_DW_AN];
// 		s_fac_new[X.k] = s_fac_sp[F_UP_CR];
// 		s_fac_new[X.k+1] = s_fac_sp[F_DW_CR];
// 		tau_new[X.k] = tau2[0]; // f+up
// 		tau_new[X.k+1] = tau2[1]; // f+down
// 		for(int i=0; i<X.k+2; i++){
// 			for(int j=i+1; j<X.k+2; j++){
// 				ln_z_fac += D0_calc_interp(K0_sp, tau_new[i], tau_new[j]) * s_fac_new[i] * s_fac_new[j];
// 			}
// 		}
// 	}
// 	double z_fac = exp( ln_z_fac );
	
	double prob = lambda_f[0] * lambda_f[1] * b_fac * z_fac * (2./ lambda_b);
	prob *= l_max_m * l_max_p * (double)Y.k / (double)((S[0].k+1) * (S[1].k+1));
	if( flag_rotate )  prob = -prob;
// 	if( UINF == 0 )  prob 
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_segsp, ACCPT_sub.n_reject_segsp) ){
		
		array_elim(Y.tau1, Y.k, j_tau);
		array_elim(Y.tau2, Y.k, j_tau);
		Y.k --;
		
		for(int s=0; s<2; s++){
			add_mat_M(S[s], S[s].k, S[s].k, Delta_tau1[s], Delta_tau2[s], lambda_f[s]);
			S[s].tau1[S[s].k] = tau1[s];
			S[s].tau2[S[s].k] = tau2[s];
			S[s].k ++;
		}
		
// 		X.op[i_tau_m] = F_UP_AN;  // S- -> f-up
// 		X.op[i_tau_p] = F_DW_AN;  // S+ -> f-down
		X.op[i_tau_m] = F_m_r;
		X.op[i_tau_p] = F_p_r;
		
// 		int i_tau_l, i_tau_r;
// 		double tau_l, tau_r;
// 		int F_l, F_r;
// 		if( i_tau_m > i_tau_p ){
// 			i_tau_l = i_tau_m+1;  i_tau_r = i_tau_p+1;
// 			tau_l = tau2[1];  tau_r = tau2[0];  // f+down f+up
// 			F_l = F_DW_CR;  F_r = F_UP_CR;
// 		}
// 		else{
// 			i_tau_l = i_tau_p+1;  i_tau_r = i_tau_m+1;
// 			tau_l = tau2[0];  tau_r = tau2[1];  // f+up f+down
// 			F_l = F_UP_CR;  F_r = F_DW_CR;
// 		}
// 		array_insert(X.tau, X.k, tau_l, i_tau_l);
// 		array_insert(X.tau, X.k+1, tau_r, i_tau_r);
// 		array_insert(X.op, X.k, F_l, i_tau_l);
// 		array_insert(X.op, X.k+1, F_r, i_tau_r);
// 		array_insert(X.state, X.k, EMP, i_tau_l);
// 		array_insert(X.state, X.k+1, EMP, i_tau_r);
		
// 		array_insert2(X.tau,   X.k, tau2[1], i_tau_m+1, tau2[0], i_tau_p+1);
		array_insert2(X.tau,   X.k, tau_m_ins, i_tau_m+1, tau_p_ins, i_tau_p+1);
		array_insert2(X.op,    X.k, F_m_l,     i_tau_m+1, F_p_l,     i_tau_p+1);
		array_insert2(X.state, X.k, state_m,   i_tau_m+1, state_p,   i_tau_p+1);
		X.k += 2;
		
		if( flag_rotate ){
			rotate_upward(X.tau, X.k);
			rotate_upward(X.op, X.k);
			rotate_upward(X.state, X.k);
		}
		
		X.tau[X.k] = X.tau[0] + prm.beta;
		X.tau[X.k+1] = X.tau[1] + prm.beta;
		X.op[X.k] = X.op[0];
		X.state[X.k] = X.state[0];
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}

static void seg2spin(int dbl_m, int dbl_p)
{
	if( S[0].k == 0 || S[1].k == 0 )  return;
	
	int F_m_l = F_DW_CR, F_m_r = F_UP_AN;  // S- -> f+dw f-up
	int F_p_l = F_UP_CR, F_p_r = F_DW_AN;  // S+ -> f+up f-dw
// 	double *tau_m_l = S[1].tau2, *tau_m_r = S[0].tau1;
// 	double *tau_p_l = S[0].tau2, *tau_p_r = S[1].tau1;
// 	if( dbl_m ){  swap(F_m_l, F_m_r);  swap(tau_m_l, tau_m_r);}
// 	if( dbl_p ){  swap(F_p_l, F_p_r);  swap(tau_p_l, tau_p_r);}
	if( dbl_m )  swap_val(F_m_l, F_m_r);
	if( dbl_p )  swap_val(F_p_l, F_p_r);
	
	//
	// choose tau_1 & tau_2
	//
// 	double tau1[2];  // f-[sigma]
// 	double tau2[2];  // f+[sigma]
	int j_tau1[2];
	int j_tau2[2];
	j_tau1[0] = rand_int(S[0].k);  // S- <- f-up
	j_tau1[1] = rand_int(S[1].k);  // S+ <- f-down
	
// 	double tau_m = S[0].tau1[j_tau1[0]];
// 	double tau_p = S[1].tau1[j_tau1[1]];
// 	double tau_m = tau_m_r[j_tau1[0]];
// 	double tau_p = tau_p_r[j_tau1[1]];
// 	int i_tau_m = tau_position(X.tau, X.k, tau_m);
// 	int i_tau_p = tau_position(X.tau, X.k, tau_p);
	int i_tau_m = tau_position(X.tau, X.k, S[0].tau1[j_tau1[0]]);
	int i_tau_p = tau_position(X.tau, X.k, S[1].tau1[j_tau1[1]]);
	
// 	if( dbl_m )  i_tau_m = (i_tau_m!=0) ? i_tau_m-1 : X.k-1;
// 	if( dbl_p )  i_tau_p = (i_tau_p!=0) ? i_tau_p-1 : X.k-1;
	if( dbl_m )  i_tau_m = (i_tau_m - 1 + X.k) % X.k;
	if( dbl_p )  i_tau_p = (i_tau_p - 1 + X.k) % X.k;
	
	#if TEST_MODE
	printf("\nseg2spin  (dbl_m=%d  dbl_p=%d)\n", dbl_m, dbl_p);
// 	printf(" j_tau1[0]=%d\n", j_tau1[0]);
// 	printf(" j_tau1[1]=%d\n", j_tau1[1]);
	printf(" i_tau_m=%d\n", i_tau_m);
	printf(" i_tau_p=%d\n", i_tau_p);
	#endif
	
// 	if( X.op[i_tau_m+1] != F_DW_CR )  return;
// 	if( X.op[i_tau_p+1] != F_UP_CR )  return;
// 	if( X.op[i_tau_m+1] != F_m_l )  return;
// 	if( X.op[i_tau_p+1] != F_p_l )  return;
	if( X.op[i_tau_m+1] != F_m_l || X.op[i_tau_m] != F_m_r )  return;
	if( X.op[i_tau_p+1] != F_p_l || X.op[i_tau_p] != F_p_r )  return;
	
	double tau_m = X.tau[i_tau_m];
	double tau_m_rem = X.tau[(i_tau_m+1)%X.k];  // f+down
	double l_max_m = X.tau[i_tau_m+2] - tau_m;
	double l_seg_m = X.tau[i_tau_m+1] - tau_m;
// 	tau2[1] = X.tau[i_tau_m+1];  // f+down
	
	double tau_p = X.tau[i_tau_p];
	double tau_p_rem = X.tau[(i_tau_p+1)%X.k];  // f+up
	double l_max_p = X.tau[i_tau_p+2] - tau_p;
	double l_seg_p = X.tau[i_tau_p+1] - tau_p;
// 	tau2[0] = X.tau[i_tau_p+1];  // f+up
	
	int flag_rotate = 0;
	if( i_tau_m+1 == X.k || i_tau_p+1 == X.k )  flag_rotate = 1;
	
// 	if( dbl_m )  swap(tau_m, tau2[1]);
// 	if( dbl_p )  swap(tau_p, tau2[0]);
	
// 	j_tau2[1] = tau_position(S[1].tau2, X.k, tau2[1]);  // f+down
// 	j_tau2[0] = tau_position(S[0].tau2, X.k, tau2[0]);  // f+up
// 	j_tau2[1] = tau_position(tau_m_l, X.k, tau_m_rem);  // f+down
// 	j_tau2[0] = tau_position(tau_p_l, X.k, tau_p_rem);  // f+up
	double tau2_1 = ( dbl_m==0 ) ? tau_m_rem : tau_m;
	double tau2_0 = ( dbl_p==0 ) ? tau_p_rem : tau_p;
	j_tau2[1] = tau_position(S[1].tau2, X.k, tau2_1);  // f+down
	j_tau2[0] = tau_position(S[0].tau2, X.k, tau2_0);  // f+up
	
// 	if( dbl_m )  swap(j_tau1[0], j_tau2[1]);
// 	if( dbl_p )  swap(j_tau1[1], j_tau2[0]);
	
	//
	// probability
	//
	double lambda_b = D0_calc_interp(D0_pm, tau_m, tau_p);
	
	// l_seg_m : DW <- EMP (if dbl_m==0),  DBL (if dbl_m==1)
	// l_seg_p : UP <- EMP (if dbl_m==0),  DBL (if dbl_m==1)
	double ef_m = ( dbl_m==0 ) ? -prm.ef[1] : prm.ef[0] + prm.U;
	double ef_p = ( dbl_p==0 ) ? -prm.ef[0] : prm.ef[1] + prm.U;
	double b_fac = exp( l_seg_m * ef_m + l_seg_p * ef_p );
// 	double b_fac = exp( -l_seg_m * prm.ef[1] - l_seg_p * prm.ef[0] );
	
	
	double lambda_f[2];
	for(int s=0; s<2; s++){
		lambda_f[s] = S[s].mat_M[j_tau1[s]][j_tau2[s]];
	}
	
	double ln_z_fac = 0;
	for(int a=0; a<2; a++){  // remove (f+up, f+down)
		if( flag_g[a] ){
	// 		double s_fac1 = s_fac_sp[F_UP_CR];  // tau2[0] f+up
	// 		double s_fac2 = s_fac_sp[F_DW_CR];  // tau2[1] f+down
			double s_fac1 = s_fac[a][F_p_l];  // tau_p_rem f+up
			double s_fac2 = s_fac[a][F_m_l];  // tau_m_rem f+down
	// 		ln_z_fac -= D0_calc_interp(K0_sp, tau2[0], tau2[1]) * s_fac1 * s_fac2;  // subtract double counting
			ln_z_fac -= D0_calc_interp(K0[a], tau_p_rem, tau_m_rem) * s_fac1 * s_fac2;  // subtract double counting
			for(int i=0; i<X.k; i++){
	// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[0]) * s_fac_sp[X.op[i]] * s_fac1;
	// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[1]) * s_fac_sp[X.op[i]] * s_fac2;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_p_rem) * s_fac[a][X.op[i]] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_m_rem) * s_fac[a][X.op[i]] * s_fac2;
			}
		}
	}
	ln_z_fac = -ln_z_fac;
	for(int a=0; a<2; a++){  // replace (f-up, f-down) with (S-, S+)
		if( flag_g[a] ){
	// 		double s_fac1 = s_fac_sp[S_MINUS] - s_fac_sp[F_UP_AN];  // tau_m : S- <- f-up
	// 		double s_fac2 = s_fac_sp[S_PLUS] - s_fac_sp[F_DW_AN];   // tau_p : S+ <- f-down
			double s_fac1 = s_fac[a][S_MINUS] - s_fac[a][F_m_r];  // tau_m : S- <- f-up
			double s_fac2 = s_fac[a][S_PLUS] - s_fac[a][F_p_r];   // tau_p : S+ <- f-down
			ln_z_fac += D0_calc_interp(K0[a], tau_p, tau_m)
			         * (s_fac[a][S_MINUS] * s_fac[a][S_PLUS] - s_fac[a][F_m_r] * s_fac[a][F_p_r]);
			
			double s_fac_cav[X.k];
			for(int i=0; i<X.k; i++)  s_fac_cav[i] = s_fac[a][X.op[i]];
			s_fac_cav[i_tau_m] = s_fac_cav[i_tau_p] = 0;
			s_fac_cav[(i_tau_m+1)%X.k] = s_fac_cav[(i_tau_p+1)%X.k] = 0;
			for(int i=0; i<X.k; i++){
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_m) * s_fac_cav[i] * s_fac1;
				ln_z_fac += D0_calc_interp(K0[a], X.tau[i], tau_p) * s_fac_cav[i] * s_fac2;
			}
		}
	}
	double z_fac = exp( ln_z_fac );
	
// 	double ln_z_fac = 0;
// 	{
// 		for(int i=0; i<X.k; i++){
// 			for(int j=i+1; j<X.k; j++){
// 				ln_z_fac -= D0_calc_interp(K0_sp, X.tau[i], X.tau[j]) * s_fac_sp[X.op[i]] * s_fac_sp[X.op[j]];
// 			}
// 		}
// 		
// 		double s_fac_new[X.k];
// 		for(int i=0; i<X.k; i++){
// 			s_fac_new[i] = s_fac_sp[X.op[i]];
// 		}
// 		s_fac_new[i_tau_m] = s_fac_sp[S_MINUS];
// 		s_fac_new[i_tau_p] = s_fac_sp[S_PLUS];
// 		s_fac_new[(i_tau_m+1)%X.k] = 0;
// 		s_fac_new[(i_tau_p+1)%X.k] = 0;
// 		for(int i=0; i<X.k; i++){
// 			for(int j=i+1; j<X.k; j++){
// 				ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], X.tau[j]) * s_fac_new[i] * s_fac_new[j];
// 			}
// 		}
// 	}
// 	double z_fac = exp( ln_z_fac );
	
	double prob = lambda_f[0] * lambda_f[1] * b_fac * z_fac * (0.5 * lambda_b);
	prob *= (double)(S[0].k * S[1].k) / (l_max_m * l_max_p * (double)(Y.k+1));
	if( flag_rotate )  prob = -prob;
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_segsp, ACCPT_sub.n_reject_segsp) ){
		
		Y.tau1[Y.k] = tau_m;
		Y.tau2[Y.k] = tau_p;
		Y.k++;
		
		for(int s=0; s<2; s++){
			remove_mat_M(S[s], j_tau1[s], j_tau2[s], lambda_f[s]);
			array_elim(S[s].tau1, S[s].k, j_tau1[s]);
			array_elim(S[s].tau2, S[s].k, j_tau2[s]);
			S[s].k --;
		}
		
		X.op[i_tau_m] = S_MINUS;  // S- <- f-up
		X.op[i_tau_p] = S_PLUS;  // S+ <- f-down
		
// 		int i_tau_l = (i_tau_m+1)%X.k;
// 		int i_tau_r = (i_tau_p+1)%X.k;
// 		if(  i_tau_l < i_tau_r ){
// 			i_tau_l = (i_tau_p+1)%X.k;
// 			i_tau_r = (i_tau_m+1)%X.k;
// 		}
		
// 		array_elim(X.tau, X.k, i_tau_l);
// 		array_elim(X.tau, X.k-1, i_tau_r);
// 		array_elim(X.op, X.k, i_tau_l);
// 		array_elim(X.op, X.k-1, i_tau_r);
// 		array_elim(X.state, X.k, i_tau_l);
// 		array_elim(X.state, X.k-1, i_tau_r);
		
// 		array_elim2(X.tau, X.k, i_tau_l, i_tau_r);
// 		array_elim2(X.op, X.k, i_tau_l, i_tau_r);
// 		array_elim2(X.state, X.k, i_tau_l, i_tau_r);
		
		array_elim2(X.tau,   X.k, (i_tau_m+1)%X.k, (i_tau_p+1)%X.k);
		array_elim2(X.op,    X.k, (i_tau_m+1)%X.k, (i_tau_p+1)%X.k);
		array_elim2(X.state, X.k, (i_tau_m+1)%X.k, (i_tau_p+1)%X.k);
		X.k -= 2;
		
		X.tau[X.k] = X.tau[0] + prm.beta;
		X.tau[X.k+1] = X.tau[1] + prm.beta;
		X.op[X.k] = X.op[0];
		X.state[X.k] = X.state[0];
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}
/*
static void seg2spin(int dbl_m, int dbl_p)
{
	if( S[0].k == 0 || S[1].k == 0 )  return;
	printf("\nseg2spin  (dbl_m=%d  dbl_p=%d)\n", dbl_m, dbl_p);
	
	//
	// choose tau_1 & tau_2
	//
	int F_m_l = F_DW_CR, F_m_r = F_UP_AN;  // S- -> f+dw f-up
	int F_p_l = F_UP_CR, F_p_r = F_DW_AN;  // S+ -> f+up f-dw
	double *tau_m_l = S[1].tau2, *tau_m_r = S[0].tau1;
	double *tau_p_l = S[0].tau2, *tau_p_r = S[1].tau1;
	if( dbl_m ){  swap(F_m_l, F_m_r);  swap(tau_m_l, tau_m_r);}
	if( dbl_p ){  swap(F_p_l, F_p_r);  swap(tau_p_l, tau_p_r);}
	
// 	double tau1[2];  // f-[sigma]
// 	double tau2[2];  // f+[sigma]
	int j_tau1[2];
	int j_tau2[2];
	j_tau1[0] = rand_int(S[0].k);  // S- <- f-up
	j_tau1[1] = rand_int(S[1].k);  // S+ <- f-down
	
// 	double tau_m = S[0].tau1[j_tau1[0]];
// 	double tau_p = S[1].tau1[j_tau1[1]];
	double tau_m = tau_m_r[j_tau1[0]];
	double tau_p = tau_p_r[j_tau1[1]];
	printf("tau_m=%.4lf  tau_p=%.4lf\n", tau_m, tau_p);
	int i_tau_m = tau_position(X.tau, X.k, tau_m);
	int i_tau_p = tau_position(X.tau, X.k, tau_p);
	
	#if TEST_MODE
	printf("\nseg2spin  (dbl_m=%d  dbl_p=%d)\n", dbl_m, dbl_p);
// 	printf(" j_tau1[0]=%d\n", j_tau1[0]);
// 	printf(" j_tau1[1]=%d\n", j_tau1[1]);
	printf(" i_tau_m=%d\n", i_tau_m);
	printf(" i_tau_p=%d\n", i_tau_p);
	#endif
	
// 	if( X.op[i_tau_m+1] != F_DW_CR )  return;
// 	if( X.op[i_tau_p+1] != F_UP_CR )  return;
	if( X.op[i_tau_m+1] != F_m_l )  return;
	if( X.op[i_tau_p+1] != F_p_l )  return;
	
	double l_max_m = X.tau[i_tau_m+2] - tau_m;
	double l_seg_m = X.tau[i_tau_m+1] - tau_m;
// 	tau2[1] = X.tau[i_tau_m+1];  // f+down
	double tau_m_rem = X.tau[i_tau_m+1];  // f+down
	
	double l_max_p = X.tau[i_tau_p+2] - tau_p;
	double l_seg_p = X.tau[i_tau_p+1] - tau_p;
// 	tau2[0] = X.tau[i_tau_p+1];  // f+up
	double tau_p_rem = X.tau[i_tau_p+1];  // f+up
	
	int flag_rotate = 0;
	if( i_tau_m+1 == X.k ){
		flag_rotate = 1;
// 		tau2[1] = X.tau[0];
		tau_m_rem = X.tau[0];
	}
	if( i_tau_p+1 == X.k ){
		flag_rotate = 1;
// 		tau2[0] = X.tau[0];
		tau_p_rem = X.tau[0];
	}
// 	j_tau2[1] = tau_position(S[1].tau2, X.k, tau2[1]);  // f+down
// 	j_tau2[0] = tau_position(S[0].tau2, X.k, tau2[0]);  // f+up
	j_tau2[1] = tau_position(tau_m_l, X.k, tau_m_rem);  // f+down
	j_tau2[0] = tau_position(tau_p_l, X.k, tau_p_rem);  // f+up
	
	if( dbl_m )  swap(j_tau1[0], j_tau2[1]);
	if( dbl_p )  swap(j_tau1[1], j_tau2[0]);
	
	//
	// probability
	//
	double lambda_b = D0_calc_interp(D0, tau_m, tau_p);
	
	// l_seg_m : DW <- EMP (if dbl_m==0),  DBL (if dbl_m==1)
	// l_seg_p : UP <- EMP (if dbl_m==0),  DBL (if dbl_m==1)
	double ef_m = ( dbl_m==0 ) ? -prm.ef[1] : prm.ef[0] + prm.U;
	double ef_p = ( dbl_p==0 ) ? -prm.ef[0] : prm.ef[1] + prm.U;
	double b_fac = exp( l_seg_m * ef_m + l_seg_p * ef_p );
// 	double b_fac = exp( -l_seg_m * prm.ef[1] - l_seg_p * prm.ef[0] );
	
	
	double lambda_f[2];
	for(int s=0; s<2; s++){
		lambda_f[s] = S[s].mat_M[j_tau1[s]][j_tau2[s]];
	}
	
	double ln_z_fac = 0;
	{  // remove (f+up, f+down)
// 		double s_fac1 = s_fac_sp[F_UP_CR];  // tau2[0] f+up
// 		double s_fac2 = s_fac_sp[F_DW_CR];  // tau2[1] f+down
		double s_fac1 = s_fac_sp[F_p_l];  // tau_p_rem f+up
		double s_fac2 = s_fac_sp[F_m_l];  // tau_m_rem f+down
// 		ln_z_fac -= D0_calc_interp(K0_sp, tau2[0], tau2[1]) * s_fac1 * s_fac2;  // subtract double counting
		ln_z_fac -= D0_calc_interp(K0_sp, tau_p_rem, tau_m_rem) * s_fac1 * s_fac2;  // subtract double counting
		for(int i=0; i<X.k; i++){
// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[0]) * s_fac_sp[X.op[i]] * s_fac1;
// 			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau2[1]) * s_fac_sp[X.op[i]] * s_fac2;
			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau_p_rem) * s_fac_sp[X.op[i]] * s_fac1;
			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau_m_rem) * s_fac_sp[X.op[i]] * s_fac2;
		}
		ln_z_fac = -ln_z_fac;
	}
	{  // replace (f-up, f-down) with (S-, S+)
// 		double s_fac1 = s_fac_sp[S_MINUS] - s_fac_sp[F_UP_AN];  // tau_m : S- <- f-up
// 		double s_fac2 = s_fac_sp[S_PLUS] - s_fac_sp[F_DW_AN];   // tau_p : S+ <- f-down
		double s_fac1 = s_fac_sp[S_MINUS] - s_fac_sp[F_m_r];  // tau_m : S- <- f-up
		double s_fac2 = s_fac_sp[S_PLUS] - s_fac_sp[F_p_r];   // tau_p : S+ <- f-down
		ln_z_fac += D0_calc_interp(K0_sp, tau_p, tau_m)
		         * (s_fac_sp[S_MINUS] * s_fac_sp[S_PLUS] - s_fac_sp[F_m_r] * s_fac_sp[F_p_r]);
		
		double s_fac_cav[X.k];
		for(int i=0; i<X.k; i++)  s_fac_cav[i] = s_fac_sp[X.op[i]];
		s_fac_cav[i_tau_m] = s_fac_cav[i_tau_p] = 0;
		s_fac_cav[(i_tau_m+1)%X.k] = s_fac_cav[(i_tau_p+1)%X.k] = 0;
		for(int i=0; i<X.k; i++){
			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau_m) * s_fac_cav[i] * s_fac1;
			ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], tau_p) * s_fac_cav[i] * s_fac2;
		}
	}
	double z_fac = exp( ln_z_fac );
	
// 	double ln_z_fac = 0;
// 	{
// 		for(int i=0; i<X.k; i++){
// 			for(int j=i+1; j<X.k; j++){
// 				ln_z_fac -= D0_calc_interp(K0_sp, X.tau[i], X.tau[j]) * s_fac_sp[X.op[i]] * s_fac_sp[X.op[j]];
// 			}
// 		}
// 		
// 		double s_fac_new[X.k];
// 		for(int i=0; i<X.k; i++){
// 			s_fac_new[i] = s_fac_sp[X.op[i]];
// 		}
// 		s_fac_new[i_tau_m] = s_fac_sp[S_MINUS];
// 		s_fac_new[i_tau_p] = s_fac_sp[S_PLUS];
// 		s_fac_new[(i_tau_m+1)%X.k] = 0;
// 		s_fac_new[(i_tau_p+1)%X.k] = 0;
// 		for(int i=0; i<X.k; i++){
// 			for(int j=i+1; j<X.k; j++){
// 				ln_z_fac += D0_calc_interp(K0_sp, X.tau[i], X.tau[j]) * s_fac_new[i] * s_fac_new[j];
// 			}
// 		}
// 	}
// 	double z_fac = exp( ln_z_fac );
	
	double prob = lambda_f[0] * lambda_f[1] * b_fac * z_fac * (0.5 * lambda_b);
	prob *= (double)(S[0].k * S[1].k) / (l_max_m * l_max_p * (double)(Y.k+1));
	if( flag_rotate )  prob = -prob;
	
	#if TEST_MODE
	if( flag_rotate )  printf(" ROTATE\n");
	printf("  prob=%.5lf\n", prob);
	if( prob < 0 )  exit(0);
	#endif
	
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_segsp, ACCPT_sub.n_reject_segsp) ){
		
		Y.tau1[Y.k] = tau_m;
		Y.tau2[Y.k] = tau_p;
		Y.k++;
		
		for(int s=0; s<2; s++){
			remove_mat_M(S[s], j_tau1[s], j_tau2[s], lambda_f[s]);
			array_elim(S[s].tau1, S[s].k, j_tau1[s]);
			array_elim(S[s].tau2, S[s].k, j_tau2[s]);
			S[s].k --;
		}
		
		X.op[i_tau_m] = S_MINUS;  // S- <- f-up
		X.op[i_tau_p] = S_PLUS;  // S+ <- f-down
		
		int i_tau_l = (i_tau_m+1)%X.k;
		int i_tau_r = (i_tau_p+1)%X.k;
		if(  i_tau_l < i_tau_r ){
			i_tau_l = (i_tau_p+1)%X.k;
			i_tau_r = (i_tau_m+1)%X.k;
		}
		
// 		array_elim(X.tau, X.k, i_tau_l);
// 		array_elim(X.tau, X.k-1, i_tau_r);
// 		array_elim(X.op, X.k, i_tau_l);
// 		array_elim(X.op, X.k-1, i_tau_r);
// 		array_elim(X.state, X.k, i_tau_l);
// 		array_elim(X.state, X.k-1, i_tau_r);
		
		array_elim2(X.tau, X.k, i_tau_l, i_tau_r);
		array_elim2(X.op, X.k, i_tau_l, i_tau_r);
		array_elim2(X.state, X.k, i_tau_l, i_tau_r);
		X.k -= 2;
		
		X.tau[X.k] = X.tau[0] + prm.beta;
		X.tau[X.k+1] = X.tau[1] + prm.beta;
		X.op[X.k] = X.op[0];
		X.state[X.k] = X.state[0];
		
		if( prob < 0 )  w_sign = -w_sign;
		
		#if TEST_MODE
		print_config();
		#endif
	}
}
*/

static void perm_boson()
{
	if( Y.k < 2 )  return;
	
// 	int j_tau_a = rand_int(Y.k);
// 	int j_tau_b = rand_int(Y.k);
// 	if(j_tau_a == j_tau_b)  return;
	
	int j_tau_a = rand_int(Y.k);
	int j_tau_b = rand_int(Y.k-1);
	if(j_tau_a == j_tau_b)  j_tau_b = Y.k-1;
	
	double tau1_a = Y.tau1[j_tau_a];  // S-
	double tau2_a = Y.tau2[j_tau_a];  // S+
	double tau1_b = Y.tau1[j_tau_b];
	double tau2_b = Y.tau2[j_tau_b];
	
	double before = D0_calc_interp(D0_pm, tau1_a, tau2_a) * D0_calc_interp(D0_pm, tau1_b, tau2_b);
	double after  = D0_calc_interp(D0_pm, tau1_a, tau2_b) * D0_calc_interp(D0_pm, tau1_b, tau2_a);
	
	double prob = after / before;
	
	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_boson, ACCPT_sub.n_reject_boson) ){
		
		Y.tau1[j_tau_a] = tau1_b;
		Y.tau1[j_tau_b] = tau1_a;
		
// 		int i_tau_a1 = tau_position(X.tau, X.k, Y.tau1[j_tau_a]);
// 		int i_tau_a2 = tau_position(X.tau, X.k, Y.tau2[j_tau_a]);
// 		int i_tau_b1 = tau_position(X.tau, X.k, Y.tau1[j_tau_b]);
// 		int i_tau_b2 = tau_position(X.tau, X.k, Y.tau2[j_tau_b]);
		
		if( prob < 0 )  w_sign = -w_sign;
	}
}


static void shift_update_X(int i_tau, double tau, int flag_rotate_up, int flag_rotate_down)
{
	X.tau[i_tau] = tau;
	
	if( flag_rotate_up ){
		rotate_upward(X.tau, X.k);
		rotate_upward(X.op, X.k);
		rotate_upward(X.state, X.k);
	}
	if( flag_rotate_down ){
		rotate_downward(X.tau, X.k);
		rotate_downward(X.op, X.k);
		rotate_downward(X.state, X.k);
	}
	
	X.tau[X.k] = X.tau[0] + prm.beta;
	X.tau[X.k+1] = X.tau[1] + prm.beta;
	X.op[X.k] = X.op[0];
	X.state[X.k] = X.state[0];
}

static void shift_tau()
{
// 	if(X.k == 0)  return;
	
	//
	// choose tau
	//
	int i_tau = rand_int(X.k);
	double tau, l_max, l_seg;
	int flag_rotate_up = 0, flag_rotate_down = 0;
	
	if( i_tau ){
		l_max = X.tau[i_tau+1] - X.tau[i_tau-1];
		tau = X.tau[i_tau-1] + rand_tau_l(l_max);
		l_seg = tau - X.tau[i_tau];
		
		if( tau >= prm.beta ){
			tau -= prm.beta;
			flag_rotate_up = 1;
		}
	}
	else{  // i_tau==0
		l_max = X.tau[1] + prm.beta - X.tau[X.k-1];
		tau = X.tau[1] - rand_tau_l(l_max);
		l_seg = tau - X.tau[0];
		if( tau < 0 ){
			tau += prm.beta;
			flag_rotate_down = 1;
		}
	}
	
	//
	// probability
	//
	double ln_z_fac = 0;
	for(int a=0; a<2; a++){
		if( flag_g[a] ){
			double s_fac1 = s_fac[a][X.op[i_tau]];
			double temp = -D0_calc_interp(K0[a], X.tau[i_tau], tau) * s_fac1;  // subtract over-sum
			for(int i=0; i<X.k; i++){
				temp += D0_calc_interp(K0[a], X.tau[i], tau) * s_fac[a][X.op[i]];  // after
				temp -= D0_calc_interp(K0[a], X.tau[i], X.tau[i_tau]) * s_fac[a][X.op[i]];  // before
			}
			ln_z_fac += s_fac1 * temp;
		}
	}
	double z_fac = exp( ln_z_fac );
	
	static int SIGMA[] = {0, 0, 1, 1, -1, -1};  // up, down, S
	static int DAG[] = {0, 1, 0, 1, 0, 1};
	
	int sigma = SIGMA[X.op[i_tau]];
	int dag = DAG[X.op[i_tau]];
	
	//
	// shift S-operator
	//
	if( sigma < 0 ){
// 		double *D_tau, *D_tau_pair;
// 		if( dag ){  // S+
// 			D_tau = Y.tau2;  D_tau_pair = Y.tau1;
// 		}
// 		else{  // S-
// 			D_tau = Y.tau1;  D_tau_pair = Y.tau2;
// 		}
// 		int j_tau = tau_position(D_tau, Y.k, X.tau[i_tau]);
// 		double lambda = D0_calc_interp(D0, D_tau_pair[j_tau], tau)  // after
// 		              / D0_calc_interp(D0, D_tau_pair[j_tau], D_tau[j_tau]);  // before
		
		int j_tau;
		double lambda;
		if( dag ){  // S+
			double tau_p = X.tau[i_tau];
			j_tau = tau_position(Y.tau2, Y.k, tau_p);
			double tau_m = Y.tau1[j_tau];
			lambda = D0_calc_interp(D0_pm, tau_m, tau)  // after
			       / D0_calc_interp(D0_pm, tau_m, tau_p);  // before
		}
		else{  // S-
			double tau_m = X.tau[i_tau];
			j_tau = tau_position(Y.tau1, Y.k, tau_m);
			double tau_p = Y.tau2[j_tau];
			lambda = D0_calc_interp(D0_pm, tau, tau_p)  // after
			       / D0_calc_interp(D0_pm, tau_m, tau_p);  // before
		}
		
		double sign_b_fac = double(1 - 2*dag); // +1, -1
		double b_fac = exp( -sign_b_fac * l_seg * (prm.ef[0] - prm.ef[1]) );
		
		double prob = lambda * b_fac * z_fac;
		
		//
		// update
		//
		if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_shift, ACCPT_sub.n_reject_shift) ){
			
			if( dag )  Y.tau2[j_tau] = tau;  // S+
			else       Y.tau1[j_tau] = tau;  // S-
			
			if( prob < 0 )  w_sign = -w_sign;
			shift_update_X(i_tau, tau, flag_rotate_up, flag_rotate_down);
		}
	}
	//
	// shift f-operator
	//
	else{
		double sign_b_fac = double(1 - 2*dag); // +1, -1
		double ef = prm.ef[sigma];
		if( occupied(X.state[i_tau], (sigma+1)%2) )  ef += prm.U;
		double b_fac = exp( -sign_b_fac * l_seg * ef );
		
		double Delta_tau[S[sigma].k];
		int n_tau;
		double lambda;
		
		if( dag ){  // f+
			for(int i=0; i<S[sigma].k; i++){
				Delta_tau[i] = G0_calc_interp(Delta[sigma], tau, S[sigma].tau1[i]);
			}
			
			n_tau = tau_position(S[sigma].tau2, S[sigma].k, X.tau[i_tau]);
			lambda = shift2_lambda(S[sigma], n_tau, Delta_tau);
		}
		else{  // f-
			for(int i=0; i<S[sigma].k; i++){
				Delta_tau[i] = G0_calc_interp(Delta[sigma], S[sigma].tau2[i], tau);
			}
			
			n_tau = tau_position(S[sigma].tau1, S[sigma].k, X.tau[i_tau]);
			lambda = shift1_lambda(S[sigma], n_tau, Delta_tau);
		}
		
		double prob = lambda * b_fac * z_fac;
		if( flag_rotate_up || flag_rotate_down )  prob = -prob;
		
		//
		// update
		//
		if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_shift, ACCPT_sub.n_reject_shift) ){
			
			if( dag ){  // f+
				S[sigma].tau2[n_tau] = tau;
				shift2_mat_M(S[sigma], n_tau, Delta_tau, lambda);
			}
			else{  // f-
				S[sigma].tau1[n_tau] = tau;
				shift1_mat_M(S[sigma], n_tau, Delta_tau, lambda);
			}
			
			if( prob < 0 )  w_sign = -w_sign;
			shift_update_X(i_tau, tau, flag_rotate_up, flag_rotate_down);
		}
	}
}

static void change_state0()
{
// 	if( X.k )  return;
	
// 	int state_before = op2state[X.op[0]];
	int state_before = X.state[0];
	int state_after;
// 	int state_after = (state_before + rand_int(2)) % 3;
	
// 	#if LIMIT_PARAM == 0
// 	#elif LIMIT_PARAM == 1
	if( flag_V ){
// 		int state_after = rand_int(3);
		if( UINF )  state_after = (rand_int(3)+1)%3;
		else  state_after = rand_int(4);
	}
// 	#elif LIMIT_PARAM == 2
	else{
		// V=0  exclude EMP(==2)
// 		int state_after = 1 - state_before;  // 0 <-> 1
// 		int state_after = rand_int(2);  // 0, 1
		// V=0  exclude EMP(==0) & DBL(==3)
// 		state_after = rand_int(2) + 1;  // 1, 2
		if( UINF )  state_after = rand_int(2) + 1;  // U = -ef = infinity
		else  state_after = rand_int(4);
	}
// 	#endif
	
	// spinless fermion (only up-spin)  exclude DW(==1)
// 		int state_after = 2 - state_before;  // 0 <-> 2
// 		int state_after = rand_int(2) * 2;  // 0, 2
	
	if( state_after == state_before )  return;
	
	// UP DW EMP
// 	double Ef[] = {prm.ef[0], prm.ef[1], 0};
	// EMP UP DW DBL
	double Ef[] = {0, prm.ef[0], prm.ef[1], prm.ef[0]+prm.ef[1]+prm.U};
	
	double prob = exp( -prm.beta * (Ef[state_after] - Ef[state_before]) );
	
	// UP DW EMP
// 	static int state2op[] = {S_MINUS, S_PLUS, F_UP_CR};
	// S- => UP
	// S+ => DW
	// f+ => empty
	
	unsigned long dummy;
	if( metropolis_abs(fabs(prob), dummy, dummy) ){
// 		X.op[0] = state2op[ state_after ];
		X.state[0] = state_after;
	}
}





//============================================================================

static void print_config()
{
	char str_op[8][6] = {"U-", "U+", "D-", "D+", "S-", "S+"};
	char str_state[8][4] = {"--", "-*", "*-", "**"};
	
	if( X.k )  printf("\n CONFIG k=%d\n", X.k);
	else       printf("\n CONFIG k=%d  %2s\n", X.k, str_state[X.state[0]]);
	
	printf("\n X.k=%d\n", X.k);
	for(int i=0; i<X.k; i++){
		printf(" %2d %s %2s %.4lf\n", i, str_op[X.op[i]], str_state[X.state[i]], X.tau[i]);
	}
	
	for(int s=0; s<2; s++){
		if( S[s].k ){
			printf("\n S[%d].k=%d\n", s, S[s].k);
			for(int i=0; i<S[s].k; i++){
				printf(" %2d %.4lf %.4lf\n", i, S[s].tau1[i], S[s].tau2[i]);
			}
		}
	}
	
	if( Y.k ){
		printf("\n Y.k=%d\n", Y.k);
		for(int i=0; i<Y.k; i++){
			printf(" %2d %.4lf %.4lf\n", i, Y.tau1[i], Y.tau2[i]);
		}
	}
}
