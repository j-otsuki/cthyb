/*

Continuous-time quantum Monte Carlo
for the impurity Anderson model (hybridization expansion)

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

static const char tag[64] = "v2.06";

#include "hyb_qmc.h"
#include "ct_qmc_share.h"
#include "matrix_update.h"
#include "operators.h"
// #include "green_func_0.h"
// #include "fft.h"

static struct hyb_qmc_params prm;
static struct num_mc n_mc;

// static struct cond_op S[N_S];
static class Operators S[N_S];
static struct green_func_0 Delta[N_S];

static struct single_particle *SP;
static struct two_particle *TP[N_S], TP_sp, TP_ch;
static struct two_particle_tr TP_tr;
static struct phys_quant PQ;
static struct phonon_g D;

static double TP_tau[N_TP+1];

static complex<double> Delta_omega[N_S][N_TAU/2];
static double moment_f[N_S];

static int w_sign;

static struct phys_quant_bin{
	long int ave_sign;
	unsigned long n_k[N_S][N_K];
	unsigned long n_ktot[N_S*N_K];
	double Gf[N_S][N_TAU+1];
	double f_number[N_S];
	int f_number_int[N_S];
	double occup_tot, occup_mom;
// 	double chi_sp[N_TP+1], chi_ch[N_TP+1];
	double chi[N_S][N_S][N_TP+1];
	double chi_tr1[2*N_TP2+1], chi_tr2[2*N_TP2+1];
	complex<double> chi_omega[N_TP_W];
} *B, *B_TOT;


// double tau[N_TAU+1];

// double delta_tau;

struct mc_accept{
	double d_accept_seg, d_reject_seg;
	double d_accept_shift, d_reject_shift;
} ACCPT;
static struct mc_accept_sub{
	unsigned long n_accept_seg, n_reject_seg;
	unsigned long n_accept_shift, n_reject_shift;
} ACCPT_sub;

static FILE *fp_log;

// If N_ADD<0 and/or N_SHIFT<0,
// N_ADD and/or N_SHIFT is optimized according to
//   N_ADD = (ave_tot_k / acceptance) * R_ADD * corr_fac
//   N_SHIFT = (2 * ave_tot_k / acceptance) * R_SHIFT * corr_fac
// The most part of the configuration is expected to be updated, when R_ADD=1 or R_SHIFT=1.
static double R_ADD;  // set by R_ADD = -N_ADD / 10
static double R_SHIFT;  // set by R_SHIFT = -N_SHIFT / 10
static const int N_ADD_MIN=N_S, N_SHIFT_MIN=1;
// correction factor: [1:MAX_R_CORR]  determined by P(k)

static int my_rank=0, process_num=1;


static void sampling(int i_measure, int n_bin, int n_sample, int n_add, int n_shift);
static void opt_n_mc(double accept_seg, double accept_shift);

static void init_mc_config();
static void init_measure();
static void init_measure_bin();

static void func_measure0();
static void func_measure1();
static void func_measure2();
static void (*func_measure[3])() = {func_measure0, func_measure1, func_measure2};

static void func_averagebin0(int);
static void func_averagebin1(int);
static void func_averagebin2(int);
static void (*func_averagebin[3])(int) = {func_averagebin0, func_averagebin1, func_averagebin2};

static void func_average0(int);
static void func_average1(int);
static void func_average2(int);
static void (*func_average[3])(int) = {func_average0, func_average1, func_average2};

static double mpi_reduce_bin(int);
static double mpi_reduce_accept();


static void add_segment(int sigma);
static void add_segment0(int sigma);
static void add_antisegment(int sigma);
static void add_antisegment0(int sigma);
static void remove_segment(int sigma);
static void remove_segment0(int sigma);
static void remove_antisegment(int sigma);
static void remove_antisegment0(int sigma);
static void state0_change(int sigma);
static void state0_change_2(int sigma1, int sigma2);
static void state0_exchange(int sigma1, int sigma2);
static void state0_change_empty_double();

// static void (* func_add_remove[4])(int)={
// 	add_segment, add_antisegment, remove_segment, remove_antisegment
// };

static void add_seg(int sigma, int anti);
static void rem_seg(int sigma, int anti);
static void (* func_add_remove[2])(int, int)={
	add_seg, rem_seg
};

static void shift_tau1(int sigma, int i_tau1);
static void shift_tau2(int sigma, int i_tau2);
static void shift_tau1_rand(int sigma);
static void shift_tau2_rand(int sigma);
static void shift_tau1_all(int sigma);
static void shift_tau2_all(int sigma);

static void (* func_shift_tau[2])(int, int)={
	shift_tau1, shift_tau2
};

static double tau_length(Operators &F);
static double tau_overlap(double, double, Operators &);



//============================================================================
//
// FUNCTIONS TO BE CALLED BY MAIN FUNCTION
//

void hybqmc_init(struct phys_quant **p_PQ, struct single_particle **p_SP, double **p_TP_tau, struct two_particle ***p_TP, struct two_particle **p_TP_sp, struct two_particle **p_TP_ch, struct two_particle_tr **p_TP_tr, struct phonon_g **p_D)
{
	// A seed of random number (determined from time if seed=0)
	unsigned long seed = RAND_SEED;

	#if HYB_QMC_MPI
	// 	MPI_Init(&argc_in, &argv_in);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &process_num);

		mt_init_by_time_mpi(&seed);
	#else
		mt_init_by_time(&seed);
	#endif // HYB_QMC_MPI

	if(my_rank==0){
		printf("\nHYBQMC_INIT  %s\n", tag);
		printf(" N_S = %d\n", N_S);
		printf(" (N_TAU, N_TP, N_TP2) = (%d, %d, %d)\n", N_TAU, N_TP, N_TP2);
		printf(" N_WARMUP = %d\n", N_WARMUP);
		printf(" PHONON = %d\n", PHONON);
		printf(" CHI_TR = %d\n", CHI_TR);
		printf(" seed = %ld\n", seed);
		printf(" MAX_R_CORR = %d\n", MAX_R_CORR);
		printf(" K_TOT_MIN = %d\n", K_TOT_MIN);
	}

	#if HYB_QMC_MPI
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
		fprintf(fp_log, "HYB-QMC INIT\n");
		fprintf(fp_log, " mt_init: seed=%ld\n", seed);
		if(HYB_QMC_MPI){
			fprintf(fp_log, " MPI parallel with %d nodes\n", process_num);
		}
		fclose(fp_log);
	}

	//
	// allocate memory
	//
	for(int s=0; s<N_S; s++)  G0_alloc(Delta[s], N_TAU);


// 	delta_tau = prm.beta / (double)N_TAU;
// 	printf("\n delta_tau = %lf\n", delta_tau);

// 	for(int i=0; i<=N_TAU; i++){
// 		tau[i] = delta_tau * (double)i;
// 	}

	SP = (single_particle *)malloc(N_S*sizeof(single_particle));

	for(int s=0; s<N_S; s++){
		TP[s] = (two_particle *)malloc(N_S*sizeof(two_particle));
	}

	B = (phys_quant_bin *)malloc(sizeof(phys_quant_bin));
	B_TOT = (phys_quant_bin *)malloc(sizeof(phys_quant_bin));


// 	*my_rank_out = my_rank;

	*p_PQ = &PQ;
	*p_SP = SP;
	*p_TP_tau = TP_tau;
	*p_TP = TP;
	*p_TP_sp = &TP_sp;
	*p_TP_ch = &TP_ch;
	*p_TP_tr = &TP_tr;
	*p_D = &D;

	//
	// default values for moment_f
	//  ave = 0
	//  curie = 1.
	//
	for(int s=0; s<N_S; s++){
		moment_f[s] = (s%2) * 2 - 1;  // -1 for even, +1 for odd
	}
	if( N_S % 2 )  moment_f[N_S-1] = 0;
}

void hybqmc_set_nmc(struct num_mc n_mc_in)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nSET_NUM_MC\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_SET_NUM_MC\n");
	}

	n_mc = n_mc_in;

	if(n_mc.N_ADD<0)  R_ADD = -(double)n_mc.N_ADD / 10.;
	else  R_ADD = 0;

	if(n_mc.N_SHIFT<0)  R_SHIFT = -(double)n_mc.N_SHIFT / 10.;
	else  R_SHIFT = 0;
}

void hybqmc_set_params(struct hyb_qmc_params prm_in)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nSET_PARAMS\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_SET_PARAMS\n");
	}

	prm = prm_in;

	#if PHONON
	for(int s1=0; s1<N_S; s1++){
		prm.ef[s1] += prm.g * prm.g / prm.w_ph;

		for(int s2=0; s2<N_S; s2++){
			prm.U[s1][s2] -= 2.0 * prm.g * prm.g / prm.w_ph;
		}
	}
	#endif // PHONON
}

//  V_sqr : integrated value of Delta(w), or iw*Delta(iw) with w->inf
void hybqmc_set_Delta(complex<double> Delta_omega_in[N_S][N_TAU/2], double V_sqr[N_S])
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nSET_DELTA\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_SET_DELTA\n");
	}

	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU/2; i++)  Delta_omega[s][i] = Delta_omega_in[s][i];

		G0_init_fft(Delta[s], Delta_omega[s], prm.beta, V_sqr[s]);
	}
}
//  Delta = V_sqr * G0
void hybqmc_set_G0(complex<double> G0_omega[N_S][N_TAU/2], double V_sqr[N_S])
{
	complex<double> (*Delta_omega_in)[N_TAU/2];  // [N_S][N_TAU/2]
	Delta_omega_in = (complex<double>(*)[N_TAU/2])malloc(sizeof(complex<double>)*N_S*N_TAU/2);
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU/2; i++)  Delta_omega_in[s][i] = G0_omega[s][i] * V_sqr[s];
	}

	hybqmc_set_Delta(Delta_omega_in, V_sqr);

	free(Delta_omega_in);
}

// [Optional]
void hybqmc_set_moment(double moment_f_in[N_S])
{
	for(int s=0; s<N_S; s++)  moment_f[s] = moment_f_in[s];

	double ave = 0, curie = 0;
	for(int s=0; s<N_S; s++){
		ave += moment_f[s];
		curie += pow(moment_f[s], 2);
	}
	ave /= double(N_S);
	curie /= double(N_S);
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_SET_MOMENT\n");
		printf(" ave = %.5lf\n", ave);
		printf(" Curie const = %.5lf\n", curie);
	}
}

void hybqmc_eval(int flag_tp)
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nEVAL_MC\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_EVAL\n");
	}
	time_t time_start = clock();

	init_mc_config();
	init_measure();

	// warming up
	if(R_SHIFT)  sampling(0, 1, N_WARMUP * process_num, 1, 1);
	else         sampling(0, 1, N_WARMUP * process_num, 1, 0);

// 	sampling(0, 1, N_WARMUP * process_num / 2, 1, 0);
// 	init_measure();
// 	if(R_SHIFT)  sampling(0, 1, N_WARMUP * process_num / 2, 1, 1);
// 	else         sampling(0, 1, N_WARMUP * process_num / 2, 1, 0);

	if( R_ADD || R_SHIFT ){
		if(my_rank==0){
			opt_n_mc(ACCPT.d_accept_seg, ACCPT.d_accept_shift);
		}
		#if HYB_QMC_MPI
			MPI_Bcast(&n_mc.N_MSR, 4, MPI_INT, 0, MPI_COMM_WORLD);
		#endif
	}

	init_measure();

	// measuring physical quantities
	sampling(flag_tp+1, n_mc.N_BIN, n_mc.N_MSR, n_mc.N_ADD, n_mc.N_SHIFT);

	if(my_rank==0 && DISPLAY){
		time_t time_end = clock();
		char str[100];
		sprint_time(str, time_end - time_start);
		printf("\n time:%s", str);
	}
}

void hybqmc_final()
{
	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nFINAL\n");
		fclose(fp_log);
	}
	if(my_rank==0 && DISPLAY){
		printf("\nHYBQMC_FINAL\n");
	}

// 	#if HYB_QMC_MPI
// 	MPI_Finalize();
// 	#endif // HYB_QMC_MPI

	for(int s=0; s<N_S; s++)  G0_free(Delta[s]);
	free(SP);
	for(int s=0; s<N_S; s++)  free(TP[s]);
	free(B);
	free(B_TOT);
}

void hybqmc_fprint_log(char *str)
{
	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "%s", str);
	fclose(fp_log);
}


//============================================================================
//
// MC SAMPLING
//

static void eval_acceptance(double n_sample, int n_add, int n_shift)
{
	double tot_seg = n_sample * (double)n_add;
	ACCPT.d_accept_seg = (double)ACCPT_sub.n_accept_seg / tot_seg;
	ACCPT.d_reject_seg = (double)ACCPT_sub.n_reject_seg / tot_seg;

	fp_log=fopen(LOG_FILE, "a");
	fprintf(fp_log, "\n segment add/rem :");
	fprintf(fp_log, " accept %.6lf, reject %.6lf\n", ACCPT.d_accept_seg, ACCPT.d_reject_seg);

	if(DISPLAY){
		printf("\n segment add/rem :");
		printf(" accept %.6lf, reject %.6lf\n", ACCPT.d_accept_seg, ACCPT.d_reject_seg);
	}

	if(n_shift){
		double tot_shift = n_sample * (double)n_shift;
		ACCPT.d_accept_shift = (double)ACCPT_sub.n_accept_shift / tot_shift;
		ACCPT.d_reject_shift = (double)ACCPT_sub.n_reject_shift / tot_shift;

		fprintf(fp_log, " segment shift   :");
		fprintf(fp_log, " accept %.6lf, reject %.6lf\n", ACCPT.d_accept_shift, ACCPT.d_reject_shift);

		if(DISPLAY){
			printf(" segment shift   :");
			printf(" accept %.6lf, reject %.6lf\n", ACCPT.d_accept_shift, ACCPT.d_reject_shift);
		}
	}
	else{
		ACCPT.d_accept_shift = 0;
		ACCPT.d_reject_shift = 0;
	}

	fclose(fp_log);
}

static void sampling(int i_measure, int n_bin, int n_sample, int n_add, int n_shift)
{
	int local_n_sample = (n_sample + process_num - 1)/(process_num);  // rounding upward
	int global_n_sample = local_n_sample * process_num;


	if(my_rank==0){
	// 	printf("\nSAMPLING\n");
	// 	printf("\nn_bin=%d  n_sample=%d  (n_add=%d  n_shift=%d)\n",
	// 	 n_bin, n_sample, n_add, n_shift);

		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\nn_bin=%d  n_sample=%d (local_n_sample=%d)  (n_add=%d  n_shift=%d)\n",
		 n_bin, global_n_sample, local_n_sample, n_add, n_shift);
		fclose(fp_log);
	}

	//
	// init acceptance
	//
	ACCPT_sub.n_accept_seg = ACCPT_sub.n_reject_seg = 0;
	ACCPT_sub.n_accept_shift = ACCPT_sub.n_reject_shift = 0;

	//
	// sampling
	//
	if(my_rank==0 && DISPLAY){
		printf("\n|---------|---------|---------|---------|---------|\n|");
	}
	int n_meter = (int)(n_bin * local_n_sample / 50);
	int i_meter = 0;

	clock_t time_start = clock();
	double time_trans_tot=0;

	for(int i=0; i<n_bin; i++){
		clock_t time_bin_start = clock();
		double time_trans_bin=0;

		init_measure_bin();

		for(int j=0; j<local_n_sample; j++){

			for(int i_add=0; i_add<n_add; i_add++){
// 				func_add_remove[rand_int(4)](rand_int(N_S));
// 				func_add_remove[rand_int(2)](rand_int(N_S), 0);
				func_add_remove[rand_int(2)](rand_int(N_S), rand_int(2));
			}

			for(int l=0; l<n_shift; l++){
				int sigma = rand_int(N_S);
				if(S[sigma].k){
					func_shift_tau[rand_int(2)](sigma, rand_int(S[sigma].k));
				}
			}

// 			for(int s=0; s<N_S; s++){
// 				if(S[s].k){
// 					for(int l=0; l<n_shift; l++){
// // 						shift_tau1(S, sigma, rand_int(S[sigma].k));
// // 						shift_tau2(S, sigma, rand_int(S[sigma].k));
// 						func_shift_tau[rand_int(2)](s, rand_int(S[s].k));
// 					}
// 				}
// 			}

			func_measure[i_measure]();

			if( ++i_meter % n_meter == 0 ){
				if(my_rank==0 && DISPLAY){
					printf("*");  fflush(stdout);
				}
			}
		}


		#if HYB_QMC_MPI
		time_trans_bin = mpi_reduce_bin(i_measure);
		time_trans_tot += time_trans_bin;
		#endif // HYB_QMC_MPI


		if(my_rank==0){
			func_averagebin[i_measure](global_n_sample);

			clock_t time_bin_end = clock();
	// 		print_time(time_start, time_end);
// 			print_time(time_bin_start, time_bin_end, LOG_FILE);
// 			print_time_mpi(time_bin_start, time_bin_end, time_trans_bin, LOG_FILE);

			char str[100];
			sprint_time_mpi(str, time_bin_start, time_bin_end, time_trans_bin);
			hybqmc_fprint_log(str);
		}
	}
	if(my_rank==0 && DISPLAY){
		printf("\n");
		fflush(stdout);
	}

	#if HYB_QMC_MPI
	time_trans_tot += mpi_reduce_accept();
	#endif // HYB_QMC_MPI

	if(my_rank==0){
		clock_t time_end = clock();
// 		print_time_mpi(time_start, time_end, time_trans_tot, LOG_FILE);

		char str[100];
		sprint_time_mpi(str, time_start, time_end, time_trans_tot);
		hybqmc_fprint_log(str);
	}

	if(my_rank==0){
		eval_acceptance(n_bin * global_n_sample, n_add, n_shift);
		func_average[i_measure](n_bin);
	}


	#if HYB_QMC_MPI
	MPI_Bcast(&ACCPT.d_accept_seg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ACCPT.d_accept_shift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

// 	printf(" %lf\n",ACCPT.d_accept_seg);
}

//
// optimize N_ADD and N_SHIFT
//
static void opt_n_mc(double accept_seg, double accept_shift)
{
	// increase # of updates when ave_k is small
// 	double ave_tot_k = PQ.ave_ktot;
	double ave_tot_k = sqrt( K_TOT_MIN * K_TOT_MIN + PQ.ave_ktot * PQ.ave_ktot );
	if( K_TOT_MIN )  ave_tot_k = MIN( ave_tot_k, (double)K_TOT_MIN * PQ.ave_ktot );
	// for high-T, to avoid a huge number of updates

	// correction factor:  [1:MAX_R_CORR]
	//  sampling number increased when P(k) is distributed around k=0
	//  and when the shape of P(k) is deviated from Gaussian
	double corr_fac = 1.;
	if( MAX_R_CORR > 1 ){
		// probability of k=0
		double p0 = 0;
		for(int s=0; s<N_S; s++)  p0 += PQ.Z_k[s][0] / double(N_S);
		double corr1 = 1.+ p0 * (double)MAX_R_CORR;  // corr1 > 1

		// deviation of P(k) from Gaussian
		double kurt = 0;
		for(int s=0; s<N_S; s++){
			double m2 = 0, m4 = 0;
			for(int i=0; i<N_K; i++){
				m2 += pow(double(i) - PQ.ave_k[s], 2) * PQ.Z_k[s][i];
				m4 += pow(double(i) - PQ.ave_k[s], 4) * PQ.Z_k[s][i];
			}
			kurt += ( m4 / (m2*m2) - 3.) / double(N_S);
		}
		double corr2 = MAX(1.+ kurt, 1.);  // corr2 > 1

		corr_fac = MIN(corr1 * corr2, (double)MAX_R_CORR);  // 1 <= corr_fac <= MAX_R_CORR
		if(my_rank==0 && DISPLAY){
			printf("\n correction factor\n");
			printf("  %.2lf * %.2lf --> %.2lf", corr1, corr2, corr_fac);
			printf("  (p0 = %.2lf, kurt = %.2lf)\n", p0, kurt);
		}
	}

	//
	// determine N_ADD if(R_ADD>0)
	//
	if( R_ADD ){
		n_mc.N_ADD = (int)(R_ADD * ave_tot_k / accept_seg * corr_fac);
		if(n_mc.N_ADD < N_ADD_MIN)  n_mc.N_ADD = N_ADD_MIN;
	}

	//
	// determine N_SHIFT if(R_SHIFT>0)
	//
	if( R_SHIFT ){
		n_mc.N_SHIFT = (int)(2.* R_SHIFT * ave_tot_k / accept_shift * corr_fac);
		if(n_mc.N_SHIFT < N_SHIFT_MIN)  n_mc.N_SHIFT = N_SHIFT_MIN;
	}

	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\n <k>_corr = %.2lf  N_ADD = %d  N_SHIFT = %d\n", ave_tot_k, n_mc.N_ADD, n_mc.N_SHIFT);
		fclose(fp_log);
		if(DISPLAY){
			printf("\n <k>_corr = %.2lf  N_ADD = %d  N_SHIFT = %d\n", ave_tot_k, n_mc.N_ADD, n_mc.N_SHIFT);
		}
	}
}


//============================================================================
//
// INITIALIZE FOR MC SAMPLING
//

static void init_mc_config()
{
// 	K=0;
	w_sign = +1;

	// for(int s=0; s<N_S; s++){
	// 	S[s].k=0;
	// 	S[s].flag=1;
	// }


// 	for(int i=0; i<=N_TP; i++)  TP_tau[i] = (double)i * prm.beta / (double)(N_TP+1);

	// [0:beta/2]
// 	for(int i=0; i<=N_TP; i++){
// 		TP_tau[i] = (double)i * prm.beta / (double)(2*N_TP);
//
// // 		printf(" %.5lf\n", TP_tau[i]);
// 	}

	double delta_tau2 = prm.beta / (double)(2*N_TP2);

	double TP2_tau[N_TP2+1];
	for(int i=0; i<=N_TP2; i++)  TP2_tau[i] = delta_tau2 * (double)i;

	// non-linear
// 	TP_tau[0] = TP2_tau[0];
// 	for(int i=1; i<=N_TP/2; i++){
// 		TP_tau[i] = TP2_tau[i];
// 		TP_tau[i+N_TP/2] = TP2_tau[i * (2*N_TP2/N_TP - 1) + N_TP/2];
// 	}

	TP_tau[0] = TP2_tau[0];
	int R_TP = N_TP2/N_TP;
	int i1=0, i2=0;
	for(int j=0; j<R_TP; j++){
		for(int i=0; i<N_TP/2/R_TP; i++){
			i1 ++;
			i2 += j+1;
			TP_tau[i1] = TP2_tau[i2];
		}
	}
	for(int j=0; j<R_TP; j++){
		for(int i=0; i<N_TP/2/R_TP; i++){
			i1 ++;
			i2 += j+R_TP;
			TP_tau[i1] = TP2_tau[i2];
		}
	}

}

static void init_measure()
{
// 	n_measure=0;

	PQ.ave_sign = 0;
	PQ.ave_sign_err = 0;

	for(int i=0; i<N_K*N_S; i++)  PQ.Z_ktot[i] = PQ.Z_ktot_err[i] = 0;
	PQ.ave_ktot = PQ.ave_ktot_err = 0;

	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_K; i++)  PQ.Z_k[s][i] = PQ.Z_k_err[s][i] = 0;
		PQ.ave_k[s] = PQ.ave_k_err[s] = 0;
	}

	//
	// SP
	//
	for(int s=0; s<N_S; s++){

		for(int i=0; i<=N_TAU; i++){
			SP[s].Gf_tau[i]=0;
// 			SP[s].Gf_tau_test[i]=0;
			SP[s].Gf_tau_err[i]=0;
// 			SP[s].n[i]=0;
		}

		SP[s].f_number = 0;
		SP[s].f_number_err = 0;
		SP[s].f_number2 = 0;
		SP[s].f_number2_err = 0;
	}

	PQ.occup_tot = PQ.occup_tot_err = 0;
	PQ.occup_mom = PQ.occup_mom_err = 0;

	//
	// TP
	//
	PQ.stat_suscep_sp=0;
	PQ.stat_suscep_ch=0;
	PQ.stat_suscep_sp_err=0;
	PQ.stat_suscep_ch_err=0;

	for(int i=0; i<=N_TP; i++){
		for(int s1=0; s1<N_S; s1++){
			for(int s2=0; s2<N_S; s2++){
				TP[s1][s2].chi_tau[i] = 0;
				TP[s1][s2].chi_tau_err[i] = 0;
			}
		}

		TP_sp.chi_tau[i] = 0;
		TP_sp.chi_tau_err[i] = 0;
		TP_ch.chi_tau[i] = 0;
		TP_ch.chi_tau_err[i] = 0;
	}

	for(int i=0; i<=2*N_TP2; i++){
		TP_tr.chi_tau1[i] = 0;
		TP_tr.chi_tau2[i] = 0;
		TP_tr.chi_tau1_err[i] = 0;
		TP_tr.chi_tau2_err[i] = 0;
	}

	for(int i=0; i<N_TP_W; i++){
		TP_tr.chi_omega[i] = 0;
		TP_tr.chi_omega_err[i] = 0;
	}
}

static void init_measure_bin()
{
	B->ave_sign = 0;

	for(int i=0; i<N_K*N_S; i++)  B->n_ktot[i] = 0;
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_K; i++)  B->n_k[s][i] = 0;
	}

	for(int s=0; s<N_S; s++){
		B->f_number[s] = 0;
		B->f_number_int[s] = 0;

		for(int i=0; i<=N_TAU; i++){
			B->Gf[s][i] = 0;
		}

		for(int s1=0; s1<N_S; s1++){
			for(int s2=0; s2<N_S; s2++){
				for(int i=0; i<=N_TP; i++)  B->chi[s1][s2][i] = 0;
			}
		}
	}
	B->occup_tot = B->occup_mom = 0;

	for(int i=0; i<=2*N_TP2; i++){
		B->chi_tr1[i] = 0;
		B->chi_tr2[i] = 0;
	}
	for(int i=0; i<N_TP_W; i++){
		B->chi_omega[i] = 0;
	}
}

//============================================================================
//
// MEASUREMENT
//

static inline void measure_stat()
{
	int K = 0;
	for(int s=0; s<N_S; s++){
		K += S[s].k;
	}
	B->n_ktot[K]++;

	for(int s=0; s<N_S; s++){
		B->n_k[s][S[s].k] ++;
	}

	B->ave_sign += w_sign;
// 	if(w_sign>0)  B->ave_sign++;
// 	else          B->ave_sign--;
}

static inline void measure_sp()
{
	//
	// f number
	//
	for(int s=0; s<N_S; s++){
// 		if(S[s].flag==0)  B->f_number_int[s] ++;
		if(S[s].flag==0)  B->f_number_int[s] += w_sign;
	}

	for(int s=0; s<N_S; s++){
		// double len = tau_length(S[s]);
		double len = S[s].length();
		B->f_number[s] += len * (double)w_sign;
		B->occup_tot += len * (double)w_sign;
		B->occup_mom += len * moment_f[s] * (double)w_sign;
	}

	//
	// Green function
	//
	double delta_tau = prm.beta / (double)N_TAU;

	for(int s=0; s<N_S; s++){
		for(int i=0; i<S[s].k; i++){
			for(int j=0; j<S[s].k; j++){
				double tau = S[s].tau1[j] - S[s].tau2[i];
				int i_tau;

				if(tau>0){
					i_tau = (int)(tau / delta_tau);
					B->Gf[s][i_tau] -= S[s].D.mat_M[j][i] * (double)w_sign;
				}
				else{
					i_tau = (int)((tau + prm.beta) / delta_tau);
					B->Gf[s][i_tau] += S[s].D.mat_M[j][i] * (double)w_sign;
				}

// 				SP[s].n[i_tau] ++;
			}
		}
	}

}


// n > 0
static inline void arrange_1array(Operators &F, double *tau, int *flag_exist, int n)
{
	// double *temp_tau1, *temp_tau2;
	std::vector<double> *temp_tau1, *temp_tau2;

	if(F.flag){
		for(int i=0; i<F.k*n; i++){
			flag_exist[2*i] = 0;
			flag_exist[2*i+1] = 1;
		}
		flag_exist[2*F.k*n] = 0;

		temp_tau1 = &F.tau1;
		temp_tau2 = &F.tau2;
	}
	else{
		for(int i=0; i<F.k*n; i++){
			flag_exist[2*i] = 1;
			flag_exist[2*i+1] = 0;
		}
		flag_exist[2*F.k*n] = 1;

		temp_tau1 = &F.tau2;
		temp_tau2 = &F.tau1;
	}

	for(int i=0; i<F.k*n; i++){
		int quot = (int)(i/F.k);
		int res = i%F.k;
		tau[2*i+1] = (*temp_tau2)[res] + prm.beta * (double)quot;
		tau[2*i+2] = (*temp_tau1)[res] + prm.beta * (double)quot;;
	}
	tau[0] = 0;
	tau[2*F.k*n+1] = prm.beta * (double)n;

	#if TEST_MODE==2
	printf("\narrange_1array  n=%d\n", n);
	for(int i=0; i<2*F.k*n+1; i++){
		printf(" %d  %8.5lf %d\n", i, tau[i], flag_exist[i]);
	}
	printf(" %d  %8.5lf\n", 2*F.k*n+1, tau[2*F.k*n+1]);
	#endif // TEST_MODE
}

static double measure_tp_sub(double TP_tau, double *temp_tau1, int *flag1, int k1, double *tau2, int *flag2, int k2)
{
	double chi_tau = 0;
// 	for(int n=0; n<=N_TP; n++){
		double tau1[2*k1+2];
		for(int i=0; i<2*k1+2; i++)  tau1[i] = temp_tau1[i] + TP_tau;

		double tau_comb[2*(k1+k2)+4];
		int flag_comb[2*(k1+k2)+4];

		int k=0, i1=0, i2=0;
		while( tau1[0] >= tau2[i2] ) i2++;

		do{
			if( tau1[i1] < tau2[i2] ){
				tau_comb[k] = tau1[i1];
				flag_comb[k] = flag1[i1] & flag2[i2-1];
				k++;
				i1++;
			}
			else{
				tau_comb[k] = tau2[i2];
				flag_comb[k] = flag1[i1-1] & flag2[i2];
				k++;
				i2++;
			}
		}while( i1 <= 2*k1+1);

		for(int i=0; i<k-1; i++){
			if(flag_comb[i]){
				chi_tau += ( tau_comb[i+1] - tau_comb[i] ) * (double)w_sign;
			}
		}

		#if TEST_MODE==2
		printf("\n tau_comb  flag_comb  TP_tau=%.5lf\n", TP_tau);
		for(int i=0; i<k; i++){
			printf(" %.5lf  %d\n", tau_comb[i], flag_comb[i]);
		}
		#endif // TEST_MODE
// 	}
	return chi_tau;
}
/*
static inline void measure_tp_sub(double *chi_tau, double *temp_tau1, int *flag1, int k1, double *tau2, int *flag2, int k2)
{

	for(int n=0; n<=N_TP; n++){
		double tau1[2*k1+2];
		for(int i=0; i<2*k1+2; i++)  tau1[i] = temp_tau1[i] + TP_tau[n];

		double tau_comb[2*(k1+k2)+4];
		int flag_comb[2*(k1+k2)+4];

		int k=0, i1=0, i2=0;
		while( tau1[0] >= tau2[i2] ) i2++;

		do{
			if( tau1[i1] < tau2[i2] ){
				tau_comb[k] = tau1[i1];
				flag_comb[k] = flag1[i1] & flag2[i2-1];
				k++;
				i1++;
			}
			else{
				tau_comb[k] = tau2[i2];
				flag_comb[k] = flag1[i1-1] & flag2[i2];
				k++;
				i2++;
			}
		}while( i1 <= 2*k1+1);

		for(int i=0; i<k-1; i++){
			if(flag_comb[i]){
				chi_tau[n] += ( tau_comb[i+1] - tau_comb[i] ) * (double)w_sign;
			}
		}

		#if TEST_MODE==2
		printf("\n tau_comb  flag_comb  TP_tau=%.5lf\n", TP_tau[n]);
		for(int i=0; i<k; i++){
			printf(" %.5lf  %d\n", tau_comb[i], flag_comb[i]);
		}
		#endif // TEST_MODE
	}
}
*/
static inline void measure_tp()
{
	double tau1[N_S][2*N_K+2], tau2[N_S][4*N_K+2];
	int flag1[N_S][2*N_K+2], flag2[N_S][4*N_K+2];

	for(int s=0; s<N_S; s++){
// 		if(S[s].k){ // !=0
			arrange_1array(S[s], tau1[s], flag1[s], 1);
			arrange_1array(S[s], tau2[s], flag2[s], 2);
// 		}
	}

	for(int s1=0; s1<N_S; s1++){
		for(int s2=0; s2<N_S; s2++){
// 			measure_tp_sub(B->chi[s1][s2], tau1[s2], flag1[s2], S[s2].k, tau2[s1], flag2[s1], S[s1].k);
			for(int n=0; n<=N_TP; n++){
				B->chi[s1][s2][n] += measure_tp_sub(TP_tau[n], tau1[s2], flag1[s2], S[s2].k, tau2[s1], flag2[s1], S[s1].k);
			}
		}
	}
}
int if_double_occup()
{
// 	void arrange_1array(Operators &F, double *tau, int *flag_exist, int n);
// 	void measure_tp_sub(double *chi_tau, double *temp_tau1, int *flag1, int k1, double *tau2, int *flag2, int k2);

	double tau1[N_S][2*N_K+2], tau2[N_S][4*N_K+2];
	int flag1[N_S][2*N_K+2], flag2[N_S][4*N_K+2];

	for(int s=0; s<N_S; s++){
		arrange_1array(S[s], tau1[s], flag1[s], 1);
		arrange_1array(S[s], tau2[s], flag2[s], 2);
	}
	{
		int s1=0, s2=1;
		double chi0 = measure_tp_sub(0, tau1[s2], flag1[s2], S[s2].k, tau2[s1], flag2[s1], S[s1].k);
		if( chi0 != 0 )  return 1;
	}
	return 0;
}

static inline void measure_tp_tr1()
{
	#if N_S==2

	double delta_tau_TP2 = prm.beta / (N_TP2*2);

	int i_tau0[S[0].k*S[1].k], m_tau0[S[0].k*S[1].k];
	int j_tau0[S[0].k*S[1].k], l_tau0[S[0].k*S[1].k];

	int i_tau1[S[0].k*S[1].k], m_tau1[S[0].k*S[1].k];
	int j_tau1[S[0].k*S[1].k], l_tau1[S[0].k*S[1].k];

	int n_ij0=0, n_lm0=0;
	int n_ij1=0, n_lm1=0;

	int sign_ij0[S[0].k*S[1].k], sign_lm0[S[0].k*S[1].k];
	int sign_ij1[S[0].k*S[1].k], sign_lm1[S[0].k*S[1].k];

	for(int i=0; i<S[0].k; i++){
		for(int j=0; j<S[1].k; j++){

			double tau = S[0].tau2[i] - S[1].tau1[j];
			int sign_ij = +1;

			if(tau<0){
				tau += prm.beta;
				sign_ij = -1;
			}
			int d_tau = (int)(tau / delta_tau_TP2);

			if(d_tau==0){
				i_tau0[n_ij0] = i;
				j_tau0[n_ij0] = j;
				sign_ij0[n_ij0++] = sign_ij;
			}
			else if(d_tau==1){
				i_tau1[n_ij1] = i;
				j_tau1[n_ij1] = j;
				sign_ij1[n_ij1++] = sign_ij;
			}
		}
	}
	for(int l=0; l<S[1].k; l++){
		for(int m=0; m<S[0].k; m++){

			double tau = S[1].tau2[l] - S[0].tau1[m];
			int sign_lm = +1;

			if(tau<0){
				tau += prm.beta;
				sign_lm = -1;
			}
			int d_tau = (int)(tau / delta_tau_TP2);

			if(d_tau==0){
				l_tau0[n_lm0] = l;
				m_tau0[n_lm0] = m;
				sign_lm0[n_lm0++] = sign_lm;
			}
			else if(d_tau==1){
				l_tau1[n_lm1] = l;
				m_tau1[n_lm1] = m;
				sign_lm1[n_lm1++] = sign_lm;
			}
		}
	}

	for(int ij=0; ij<n_ij0; ij++){
		for(int lm=0; lm<n_lm0; lm++){
			double tau = 0.5 * (S[1].tau1[j_tau0[ij]] - S[0].tau1[m_tau0[lm]]
			                  + S[0].tau2[i_tau0[ij]] - S[1].tau2[l_tau0[lm]]);

			if(tau<0){
				tau += prm.beta;
			}
			int d_tau = (int)(tau / delta_tau_TP2);

// 			B->chi_tr0[d_tau] -= S[1].D.mat_M[j_tau0[ij]][l_tau0[lm]] * S[0].D.mat_M[m_tau0[lm]][i_tau0[ij]]
// 			                   * (double)(sign_ij0[ij] * sign_lm0[lm]);

			double temp_chi_tr = -S[1].D.mat_M[j_tau0[ij]][l_tau0[lm]] * S[0].D.mat_M[m_tau0[lm]][i_tau0[ij]]
			                    * (double)(sign_ij0[ij] * sign_lm0[lm]);

			B->chi_tr1[d_tau] += 1.5 * temp_chi_tr * (double)w_sign;
		}
	}
	for(int ij=0; ij<n_ij1; ij++){
		for(int lm=0; lm<n_lm1; lm++){
			double tau = 0.5 * (S[1].tau1[j_tau1[ij]] - S[0].tau1[m_tau1[lm]]
			                  + S[0].tau2[i_tau1[ij]] - S[1].tau2[l_tau1[lm]]);

			if(tau<0){
				tau += prm.beta;
			}
			int d_tau = (int)(tau / delta_tau_TP2);

// 			B->chi_tr1[d_tau] -= S[1].D.mat_M[j_tau1[ij]][l_tau1[lm]] * S[0].D.mat_M[m_tau1[lm]][i_tau1[ij]]
// 			                   * (double)(sign_ij1[ij] * sign_lm1[lm]);

			double temp_chi_tr = -S[1].D.mat_M[j_tau1[ij]][l_tau1[lm]] * S[0].D.mat_M[m_tau1[lm]][i_tau1[ij]]
			                    * (double)(sign_ij1[ij] * sign_lm1[lm]);

			B->chi_tr1[d_tau] += -0.5 * temp_chi_tr * (double)w_sign;
		}
	}

	#endif
}


static inline void measure_tp_tr2()
{
	#if N_S==2

	double delta_tau_TP2 = prm.beta / (N_TP2*2);

	int n_ij=0, n_lm=0;

	int i_tau[S[0].k*S[1].k], m_tau[S[0].k*S[1].k];
	int j_tau[S[0].k*S[1].k], l_tau[S[0].k*S[1].k];

	double tau_ij[S[0].k*S[1].k], tau_lm[S[0].k*S[1].k];

	int sign_ij[S[0].k*S[1].k], sign_lm[S[0].k*S[1].k];

	for(int i=0; i<S[0].k; i++){
		for(int j=0; j<S[1].k; j++){

			double tau = S[0].tau2[i] - S[1].tau1[j];
			int sign = +1;

			if(tau<0){
				tau += prm.beta;
				sign = -1;
			}

			if( (int)(tau / delta_tau_TP2) == 0 ){
				i_tau[n_ij] = i;
				j_tau[n_ij] = j;
				tau_ij[n_ij] = tau / delta_tau_TP2;
				sign_ij[n_ij++] = sign;
			}
		}
	}
	for(int l=0; l<S[1].k; l++){
		for(int m=0; m<S[0].k; m++){

			double tau = S[1].tau2[l] - S[0].tau1[m];
			int sign = +1;

			if(tau<0){
				tau += prm.beta;
				sign = -1;
			}

			if( (int)(tau / delta_tau_TP2) == 0 ){
				l_tau[n_lm] = l;
				m_tau[n_lm] = m;
				tau_lm[n_lm] = tau / delta_tau_TP2;
				sign_lm[n_lm++] = sign;
			}
		}
	}

	double two_pi_t = 2.0 * M_PI / prm.beta;

	for(int ij=0; ij<n_ij; ij++){
		for(int lm=0; lm<n_lm; lm++){
			double tau = 0.5 * (S[1].tau1[j_tau[ij]] - S[0].tau1[m_tau[lm]]
			                  + S[0].tau2[i_tau[ij]] - S[1].tau2[l_tau[lm]]);

			if(tau<0){
				tau += prm.beta;
			}
			int d_tau = (int)(tau / delta_tau_TP2);

// 			B->chi_tr[d_tau] -= S[1].D.mat_M[j_tau[ij]][l_tau[lm]] * S[0].D.mat_M[m_tau[lm]][i_tau[ij]]
// 			                   * (double)(sign_ij[ij] * sign_lm[lm]);

			double temp_chi_tr = -S[1].D.mat_M[j_tau[ij]][l_tau[lm]] * S[0].D.mat_M[m_tau[lm]][i_tau[ij]]
			                   * (double)(sign_ij[ij] * sign_lm[lm]);


			double weight1 = 0;
			if( tau_ij[ij] < 0.5 ){
				if( tau_lm[lm] < 0.5 )  weight1 = 6.0;
			}
			else{
				if( tau_lm[lm] > 0.5 ) weight1 = -2.0;
			}
			B->chi_tr1[d_tau] += weight1 * temp_chi_tr * (double)w_sign;


			double weight2 = -6.0 * (tau_ij[ij] + tau_lm[lm]) + 7.0;
			complex<double> temp = weight2 * temp_chi_tr;
			B->chi_tr2[d_tau] += real(temp) * (double)w_sign;


			complex<double> fac = exp(IMAG * two_pi_t * tau);
			B->chi_omega[0] += temp * (double)w_sign;
			for(int i=1; i<N_TP_W; i++){
				temp *= fac;
				B->chi_omega[i] += temp * (double)w_sign;
			}
		}
	}

	#endif
}


static void func_measure0()
{
	measure_stat();
}

static void func_measure1()
{
	measure_stat();
	measure_sp();
}

static void func_measure2()
{
	measure_stat();
	measure_sp();
	measure_tp();

	#if CHI_TR
// 	measure_tp_tr1();
	measure_tp_tr2();
	#endif // CHI_TR
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


	unsigned long sum_n_k = 0;
	for(int i=0; i<N_K; i++)  sum_n_k += B->n_k[0][i];

	double B_ave_ktot = 0;
	for(int i=0; i<N_S*N_K; i++){
		double B_n_ktot = double(B->n_ktot[i]) / double(sum_n_k);
		PQ.Z_ktot[i] += B_n_ktot;
		PQ.Z_ktot_err[i] += pow(B_n_ktot, 2);

		B_ave_ktot += B_n_ktot * double(i);
	}
	PQ.ave_ktot += B_ave_ktot;
	PQ.ave_ktot_err += pow(B_ave_ktot, 2);

	double B_ave_k[N_S] = {0};
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_K; i++){
			double B_n_k = double(B->n_k[s][i]) / double(sum_n_k);
			PQ.Z_k[s][i] += B_n_k;
			PQ.Z_k_err[s][i] += pow(B_n_k, 2);

			B_ave_k[s] += B_n_k * double(i);
		}
		PQ.ave_k[s] += B_ave_k[s];
		PQ.ave_k_err[s] += pow(B_ave_k[s], 2);
	}
}

static void averagebin_sp(int n_sample)
{

	for(int s=0; s<N_S; s++){
		B->f_number[s] /= prm.beta * (double)n_sample;

		SP[s].f_number += B->f_number[s];
		SP[s].f_number_err += pow(B->f_number[s], 2);
	}

	for(int s=0; s<N_S; s++){
		double B_f_number2 = (double)B->f_number_int[s] / (double)n_sample;

		SP[s].f_number2 += B_f_number2;
		SP[s].f_number2_err += pow(B_f_number2, 2);
	}

	B->occup_tot /= prm.beta * (double)n_sample;
	B->occup_mom /= prm.beta * (double)n_sample;
	PQ.occup_tot += B->occup_tot;
	PQ.occup_mom += B->occup_mom;
	PQ.occup_tot_err += pow(B->occup_tot, 2);
	PQ.occup_mom_err += pow(B->occup_mom, 2);


	double delta_tau = prm.beta / (double)N_TAU;
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU; i++){
			B->Gf[s][i] /= prm.beta * delta_tau * (double)n_sample;
		}

		for(int i=N_TAU-1; i>0; i--){
			B->Gf[s][i] = (B->Gf[s][i] + B->Gf[s][i-1]) * 0.5;
		}

// 		B->Gf[s][0] = B->f_number[s] - 1.0;
// 		B->Gf[s][N_TAU] = - B->f_number[s];

// 		for(int i=0; i<=N_TAU; i++){
		for(int i=1; i<N_TAU; i++){
			SP[s].Gf_tau[i] += B->Gf[s][i];
			SP[s].Gf_tau_err[i] += pow(B->Gf[s][i], 2);
		}
	}

}

static void averagebin_tp(int n_sample)
{
	for(int s1=0; s1<N_S; s1++){
		for(int s2=0; s2<N_S; s2++){
			for(int i=0; i<=N_TP; i++){
				B->chi[s1][s2][i] /= prm.beta * (double)n_sample;
			}

			if(s1==s2)  B->chi[s1][s1][0] = B->f_number[s1];

			double temp = B->f_number[s1] * B->f_number[s2];
			for(int i=0; i<=N_TP; i++)  B->chi[s1][s2][i] -= temp;

			for(int i=0; i<=N_TP; i++){
				TP[s1][s2].chi_tau[i] += B->chi[s1][s2][i];
				TP[s1][s2].chi_tau_err[i] += pow(B->chi[s1][s2][i], 2);
			}
		}
	}

	double fac[N_S][N_S];
	for(int s1=0; s1<N_S; s1++){
		for(int s2=0; s2<N_S; s2++)  fac[s1][s2] = moment_f[s1] * moment_f[s2];
	}

	double B_chi_sp[N_TP+1] = {0}, B_chi_ch[N_TP+1];
	for(int i=0; i<=N_TP; i++){
		double B_chi_diag = 0, B_chi_offd = 0;
		for(int s1=0; s1<N_S; s1++){
			for(int s2=0; s2<N_S; s2++){
				if(s1==s2)  B_chi_diag += B->chi[s1][s2][i];
				else        B_chi_offd += B->chi[s1][s2][i];
			}
		}
		B_chi_ch[i] = B_chi_diag + B_chi_offd;

// 		B_chi_sp[i] = B_chi_diag - B_chi_offd / (double)(N_S-1);
		for(int s1=0; s1<N_S; s1++){
			for(int s2=0; s2<N_S; s2++)  B_chi_sp[i] += fac[s1][s2] * B->chi[s1][s2][i];
		}

		TP_sp.chi_tau[i] += B_chi_sp[i];
		TP_ch.chi_tau[i] += B_chi_ch[i];

		TP_sp.chi_tau_err[i] += pow(B_chi_sp[i], 2);
		TP_ch.chi_tau_err[i] += pow(B_chi_ch[i], 2);
	}


	//
	// static magnetic susceptibility
	//

	// trapezoidal rule
// 	double B_stat_suscep_sp = 0.5 * (B_chi_sp[0] + B_chi_sp[N_TP]);
// 	double B_stat_suscep_ch = 0.5 * (B_chi_ch[0] + B_chi_ch[N_TP]);
// 	for(int i=1; i<N_TP; i++){
// 		B_stat_suscep_sp += B_chi_sp[i];
// 		B_stat_suscep_ch += B_chi_ch[i];
// 	}
// 	B_stat_suscep_sp *= prm.beta / (double)(2*N_TP) * 2.0;  // extend to [0:beta]
// 	B_stat_suscep_ch *= prm.beta / (double)(2*N_TP) * 2.0;  // extend to [0:beta]

	// trapezoidal rule (non-linear mesh)
	double B_stat_suscep_sp = 0;
	double B_stat_suscep_ch = 0;
	for(int i=0; i<N_TP; i++){
		B_stat_suscep_sp += (B_chi_sp[i+1] + B_chi_sp[i]) * (TP_tau[i+1] - TP_tau[i]);
		B_stat_suscep_ch += (B_chi_ch[i+1] + B_chi_ch[i]) * (TP_tau[i+1] - TP_tau[i]);
	}
	// extend to [0:beta] (factor 0.5 was omitted instead)

	PQ.stat_suscep_sp += B_stat_suscep_sp;
	PQ.stat_suscep_ch += B_stat_suscep_ch;
	PQ.stat_suscep_sp_err += pow(B_stat_suscep_sp, 2);
	PQ.stat_suscep_ch_err += pow(B_stat_suscep_ch, 2);


	//
	// transverse susceptibility
	//
	double delta_tau_TP2 = prm.beta / (N_TP2*2);

	for(int i=0; i<2*N_TP2; i++){
		B->chi_tr1[i] /= prm.beta * pow(delta_tau_TP2,3) * (double)n_sample;
		B->chi_tr2[i] /= prm.beta * pow(delta_tau_TP2,3) * (double)n_sample;
	}
	for(int i=2*N_TP2-1; i>0; i--){
		B->chi_tr1[i] = (B->chi_tr1[i] + B->chi_tr1[i-1]) * 0.5;
		B->chi_tr2[i] = (B->chi_tr2[i] + B->chi_tr2[i-1]) * 0.5;
	}
	for(int i=0; i<=2*N_TP2; i++){
		TP_tr.chi_tau1[i] += B->chi_tr1[i];
		TP_tr.chi_tau2[i] += B->chi_tr2[i];

		TP_tr.chi_tau1_err[i] += pow(B->chi_tr1[i], 2);
		TP_tr.chi_tau2_err[i] += pow(B->chi_tr2[i], 2);
	}

	TP_tr.chi_tau1[0] = 0;
	TP_tr.chi_tau1[2*N_TP2] = 0;
	TP_tr.chi_tau1_err[0] = 0;
	TP_tr.chi_tau1_err[2*N_TP2] = 0;

	TP_tr.chi_tau2[0] = 0;
	TP_tr.chi_tau2[2*N_TP2] = 0;
	TP_tr.chi_tau2_err[0] = 0;
	TP_tr.chi_tau2_err[2*N_TP2] = 0;

	for(int i=0; i<N_TP_W; i++){
		B->chi_omega[i] /= prm.beta * pow(delta_tau_TP2,2) * (double)n_sample;

		TP_tr.chi_omega[i] += B->chi_omega[i];
		TP_tr.chi_omega_err[i] += pow(real(B->chi_omega[i]), 2) + IMAG * pow(imag(B->chi_omega[i]), 2);
	}
}

static void func_averagebin0(int n_sample)
{
	averagebin_stat(n_sample);
}
static void func_averagebin1(int n_sample)
{
	averagebin_stat(n_sample);
	averagebin_sp(n_sample);
}

static void func_averagebin2(int n_sample)
{
	averagebin_stat(n_sample);
	averagebin_sp(n_sample);
	averagebin_tp(n_sample);
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
// static inline void average_sub(double &msr, double &msr_err, int n_bin)
// {
// 	msr /= (double)n_bin;
// 	msr_err /= (double)n_bin;
// 	msr_err = sqrt(msr_err - msr*msr) * sqrt( (double)n_bin / (double)(n_bin-1) );
//
// }
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
// static inline void average_sub(complex<double> &msr, complex<double> &msr_err, int n_bin)
// {
// 	msr /= (double)n_bin;
// 	msr_err /= (double)n_bin;
//
// 	double fac_err = sqrt( (double)n_bin / (double)(n_bin-1) );
// 	msr_err = sqrt(real(msr_err) - pow(real(msr), 2)) * fac_err
// 	        + sqrt(imag(msr_err) - pow(imag(msr), 2)) * fac_err * IMAG;
// }

static void fft_chi_after_interp(double *chi_tau, complex<double> *chi_omega)
{
	double chi_tau2[2*N_TP2+1];

	// interpolation
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (INTERP, N_TP+1);
	gsl_spline_init(spline, TP_tau, chi_tau, N_TP+1);

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
// static void fft_chi_after_double(double *chi_tau, complex<double> *chi_omega)
// {
// 	double chi_tau_double[2*N_TP+1];
//
// 	// extend range from [0:beta/2] to [0:beta]
// 	for(int i=0; i<=N_TP; i++){
//
// 		chi_tau_double[i] = chi_tau[i];
// 		chi_tau_double[2*N_TP-i] = chi_tau[i];
// 	}
// 	fft_boson_radix2_tau2omega(chi_tau_double, chi_omega, prm.beta, 2*N_TP);
// }

static void average_stat(int n_bin)
{
// 	average_sub(PQ.ave_sign, PQ.ave_sign_err, n_bin);
	average_sub(PQ.ave_sign, PQ.ave_sign_err, n_bin, 1.0);

	if(my_rank==0 && DISPLAY){
		printf("\n average sign\n");
		printf("   %.6lf +- %.6lf", PQ.ave_sign, PQ.ave_sign_err);
		printf("  (deviation from 1:");
		printf(" %.3e +- %.3e)\n", 1.0-PQ.ave_sign, PQ.ave_sign_err);
	}


	for(int i=0; i<N_K*N_S; i++){
		average_sub(PQ.Z_ktot[i], PQ.Z_ktot_err[i], n_bin, PQ.ave_sign);
	}
	average_sub(PQ.ave_ktot, PQ.ave_ktot_err, n_bin, PQ.ave_sign);

	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_K; i++){
			average_sub(PQ.Z_k[s][i], PQ.Z_k_err[s][i], n_bin, PQ.ave_sign);
		}
		average_sub(PQ.ave_k[s], PQ.ave_k_err[s], n_bin, PQ.ave_sign);
	}

	int max_tot_nk;
	for(int i=1; i<N_S*N_K; i++){
		if( B->n_ktot[i] != 0 )  max_tot_nk = i;
	}
	int max_nk[N_S];
	for(int s=0; s<N_S; s++){
		for(int i=1; i<N_K; i++){
			if( PQ.Z_k[s][i] != 0 )  max_nk[s] = i;
		}
	}

	if(my_rank==0){
		fp_log=fopen(LOG_FILE, "a");
		fprintf(fp_log, "\n distribution of k");
		fprintf(fp_log, "\n  ave. :  %5.1lf :", PQ.ave_ktot);
		for(int s=0; s<N_S; s++)  fprintf(fp_log, " %5.1lf", PQ.ave_k[s]);
// 		fprintf(fp_log, "\n  peak : %3d   :", peak_tot_nk);
// 		for(int s=0; s<N_S; s++)  fprintf(fp_log, " %3d  ", peak_nk[s]);
		fprintf(fp_log, "\n  max  : %4d   :", max_tot_nk);
		for(int s=0; s<N_S; s++)  fprintf(fp_log, " %3d  ", max_nk[s]);
		fprintf(fp_log, "\n");
		fclose(fp_log);

		if(DISPLAY){
			printf("\n distribution of k");
			printf("\n  ave. :  %5.1lf :", PQ.ave_ktot);
			for(int s=0; s<N_S; s++)  printf(" %5.1lf", PQ.ave_k[s]);
// 			printf("\n  peak : %3d   :", peak_tot_nk);
// 			for(int s=0; s<N_S; s++)  printf(" %3d  ", peak_nk[s]);
			printf("\n  max  : %4d   :", max_tot_nk);
			for(int s=0; s<N_S; s++)  printf(" %3d  ", max_nk[s]);
			printf("\n");
		}
	}
}

static void average_sp(int n_bin)
{
	for(int s=0; s<N_S; s++){
		average_sub(SP[s].f_number, SP[s].f_number_err, n_bin, PQ.ave_sign);
		average_sub(SP[s].f_number2, SP[s].f_number2_err, n_bin, PQ.ave_sign);
	}
	average_sub(PQ.occup_tot, PQ.occup_tot_err, n_bin, PQ.ave_sign);
	average_sub(PQ.occup_mom, PQ.occup_mom_err, n_bin, PQ.ave_sign);

	for(int s=0; s<N_S; s++){
// 		for(int i=0; i<=N_TAU; i++){
		for(int i=1; i<N_TAU; i++){
			average_sub(SP[s].Gf_tau[i], SP[s].Gf_tau_err[i], n_bin, PQ.ave_sign);
		}

		SP[s].Gf_tau[N_TAU] = -SP[s].f_number;
		SP[s].Gf_tau_err[N_TAU] = SP[s].f_number_err;

		if( prm.UINF ){
			SP[s].Gf_tau[0] = PQ.occup_tot - 1.0;
			SP[s].Gf_tau_err[0] = PQ.occup_tot_err;
			SP[s].jump = - SP[s].Gf_tau[N_TAU] - SP[s].Gf_tau[0];
		}
		else{
			SP[s].Gf_tau[0] = SP[s].f_number - 1.0;
			SP[s].Gf_tau_err[0] = SP[s].f_number_err;
			SP[s].jump = 1;
		}
// 		printf("## jump=%lf\n", SP[s].jump);
	}

	//
	// fft
	//
	for(int s=0; s<N_S; s++){
		fft_fermion_radix2_tau2omega(SP[s].Gf_tau, SP[s].Gf_omega, prm.beta, N_TAU, SP[s].jump);
	}

	// self-energy
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU/2; i++){
			complex<double> i_omega_f = IMAG * (double)(2*i+1) * M_PI / prm.beta;

// 			SP[s].Gf0_omega[i] = 1.0 / (i_omega_f - prm.ef[s] - Delta_omega[i]);
// 			SP[s].Gf0_omega[i] = 1.0 / (i_omega_f - prm.ef[s] - Delta_omega[i] - prm.U*0.5);

			SP[s].self_f[i] = i_omega_f - prm.ef[s] - Delta_omega[s][i] - 1.0 / SP[s].Gf_omega[i];
// 			SP[s].self_f[i] = i_omega_f - prm.ef[s] - Delta_omega[i] - prm.U*0.5 - 1.0 / SP[s].Gf_omega[i];
// 			SP[s].self_f[i] = ((i_omega_f - prm.ef[s] - Delta_omega[i] - prm.U*0.5)*SP[s].Gf_omega[i] - 1.0) / SP[s].Gf_omega[i];

			#if PHONON
			SP[s].self_f[i] += prm.g * prm.g / prm.w_ph;
			#endif // PHONON

		}
	}

}

static void average_tp(int n_bin)
{

	for(int s1=0; s1<N_S; s1++){
		for(int s2=0; s2<N_S; s2++){
			for(int i=0; i<=N_TP; i++){
				average_sub(TP[s1][s2].chi_tau[i], TP[s1][s2].chi_tau_err[i], n_bin, PQ.ave_sign);
			}
		}
	}

	for(int i=0; i<=N_TP; i++){
		average_sub(TP_sp.chi_tau[i], TP_sp.chi_tau_err[i], n_bin, PQ.ave_sign);
		average_sub(TP_ch.chi_tau[i], TP_ch.chi_tau_err[i], n_bin, PQ.ave_sign);
	}

	average_sub(PQ.stat_suscep_sp, PQ.stat_suscep_sp_err, n_bin, PQ.ave_sign);
	average_sub(PQ.stat_suscep_ch, PQ.stat_suscep_ch_err, n_bin, PQ.ave_sign);


	//
	// fft
	//
	fft_chi_after_interp(TP_sp.chi_tau, TP_sp.chi_omega);
	fft_chi_after_interp(TP_ch.chi_tau, TP_ch.chi_omega);


	// checking symmetry
	for(int i=0; i<N_TP2; i++){
		if( fabs(imag(TP_sp.chi_omega[i])) > 1.0e-10 ||
		    fabs(imag(TP_ch.chi_omega[i])) > 1.0e-10 ){
			printf("\n*** error: symmetry of chi_sp or chi_ch\n");
			exit(0);
		}
	}

	// checking static susceptibility
	if( fabs(real(TP_sp.chi_omega[0]) - PQ.stat_suscep_sp) > 1.0e-10 ||
	    fabs(real(TP_ch.chi_omega[0]) - PQ.stat_suscep_ch) > 1.0e-10 ){
		printf("\n*** error: static susceptibility\n");
		exit(0);
	}


	for(int i=1; i<2*N_TP2; i++){
		average_sub(TP_tr.chi_tau1[i], TP_tr.chi_tau1_err[i], n_bin, PQ.ave_sign);
		average_sub(TP_tr.chi_tau2[i], TP_tr.chi_tau2_err[i], n_bin, PQ.ave_sign);
	}

	fft_boson_radix2_tau2omega(TP_tr.chi_tau1, TP_tr.chi_omega1, prm.beta, 2*N_TP2);
	fft_boson_radix2_tau2omega(TP_tr.chi_tau2, TP_tr.chi_omega2, prm.beta, 2*N_TP2);


	for(int i=0; i<N_TP_W; i++){
		average_sub(TP_tr.chi_omega[i], TP_tr.chi_omega_err[i], n_bin, PQ.ave_sign);
	}


	//
	// phonon Green function
	//
	#if PHONON

	double nf = 0;
	for(int s=0; s<N_S; s++)  nf += SP[s].f_number;

	{  // i=0
		double d0 = -2.0 / prm.w_ph;

		complex<double> temp = prm.g*prm.g * d0*d0 * (TP_ch.chi_omega[0] + (nf-1.0)*(nf-1.0));

		D.d_omega[0] = d0 - temp;
		D.occup = real(temp) * 0.5;
	}

	for(int i=1; i<N_TP2; i++){  // i>0
		double omega_n = (double)(2*i) * M_PI / prm.beta;
		double d0 = -2.0 * prm.w_ph / (omega_n * omega_n + prm.w_ph * prm.w_ph);

		complex<double> temp = prm.g*prm.g * d0*d0 * (TP_ch.chi_omega[i] + (nf-1.0)*(nf-1.0));

		D.d_omega[i] = d0 - temp;
		D.occup += real(temp);
	}

	D.occup /= prm.beta;
	D.occup += 1.0 / (exp(prm.beta*prm.w_ph) - 1.0);

	#endif // PHONON

}

static void func_average0(int n_bin)
{
	average_stat(n_bin);
}
static void func_average1(int n_bin)
{
	average_stat(n_bin);
	average_sp(n_bin);
}
static void func_average2(int n_bin)
{
	average_stat(n_bin);
	average_sp(n_bin);
	average_tp(n_bin);
}
//============================================================================
//
// MPI REDUCE
//

// return time
static double mpi_reduce_bin(int i_measure)
{
	double start=0, end=0;  // measure time

	#if HYB_QMC_MPI

// 	static struct phys_quant_bin{
// 		long int ave_sign;
// 		unsigned long n_k[N_S][N_K];
// 		unsigned long n_ktot[N_S*N_K];
// 		double Gf[N_S][N_TAU+1];
// 		double f_number[N_S];
// 		int f_number_int[N_S];
// 	// 	double chi_sp[N_TP+1], chi_ch[N_TP+1];
// 		double chi[N_S][N_S][N_TP+1];
// 		double chi_tr1[2*N_TP2+1], chi_tr2[2*N_TP2+1];
// 		complex<double> chi_omega[N_TP_W];
// 	} *B, *B_TOT;

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	MPI_Reduce(&B->ave_sign, &B_TOT->ave_sign, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&B->n_k[0][0], &B_TOT->n_k[0][0], N_S*N_K, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&B->n_ktot[0], &B_TOT->n_ktot[0], N_S*N_K, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if(i_measure>0){

		MPI_Reduce(&B->Gf[0][0], &B_TOT->Gf[0][0], N_S*(N_TAU+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->f_number[0], &B_TOT->f_number[0], N_S, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->f_number_int[0], &B_TOT->f_number_int[0], N_S, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->occup_tot, &B_TOT->occup_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->occup_mom, &B_TOT->occup_mom, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	if( i_measure>1 ){
// 		MPI_Reduce(&B->chi_sp[0], &B_TOT->chi_sp[0], N_TP+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		MPI_Reduce(&B->chi_ch[0], &B_TOT->chi_ch[0], N_TP+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->chi[0][0][0], &B_TOT->chi[0][0][0], N_S*N_S*(N_TP+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		#if CHI_TR
		MPI_Reduce(&B->chi_tr1[0], &B_TOT->chi_tr1[0], (2*N_TP2+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&B->chi_tr2[0], &B_TOT->chi_tr2[0], (2*N_TP2+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce(&B->chi_omega[0], &B_TOT->chi_omega[0], (2*N_TP_W), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		#endif // CHI_TR
	}

	if(my_rank==0)  *B = *B_TOT;

	end = MPI_Wtime();

// 	if(my_rank==0)  print_time_mpi(end - start);

	#endif // HYB_QMC_MPI

	return(end - start);

}

static double mpi_reduce_accept()
{
	double start=0, end=0;  // measure time

	#if HYB_QMC_MPI

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();


// 	unsigned long n_accept_seg, n_reject_seg;
// 	unsigned long n_accept_shift, n_reject_shift;

	struct mc_accept_sub ACCPT_sub_tot;
	MPI_Reduce(&ACCPT_sub.n_accept_seg, &ACCPT_sub_tot.n_accept_seg, 4, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank==0){
		ACCPT_sub = ACCPT_sub_tot;
	}

	end = MPI_Wtime();

// 	if(my_rank==0)  print_time_mpi(end - start);

	#endif // HYB_QMC_MPI

	return(end - start);
}


//============================================================================
//
// UPDATES OF MC CONFIGURATIONS
//
//============================================================================


static void add_segment0(int sigma)
{
	Operators *F = &S[sigma];

	if(F->flag==0)  return;  // |1>

// 	Operators *F_ref[N_S-1];
	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
// 		F_ref[r] = &S[s_ref[r]];
	}


	F->tau1[0] = rand_tau_beta(prm.beta);
	F->tau2[0] = rand_tau_beta(prm.beta);

	#if TEST_MODE
	printf("\nadd_segment0  (sigma=%d)\n", sigma);
	printf("  0  %8.5lf  %8.5lf\n", F->tau1[0], F->tau2[0]);
	#endif

	double lambda = G0_calc_interp(Delta[sigma], F->tau2[0], F->tau1[0]);


	double l_seg = F->tau1[0] - F->tau2[0];
	if(l_seg<0)  l_seg += prm.beta;

	double ln_b_fac = -l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0
			double l_over = tau_overlap(F->tau2[0], F->tau1[0], S[s_ref[r]]);
			ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac -= l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(F->tau1[0] - F->tau2[0]))
		                 - exp(-prm.w_ph * (prm.beta - fabs(F->tau1[0] - F->tau2[0])));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau1[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau1[0])))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau2[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau2[0])))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau1[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau1[0])))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau2[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau2[0])));
			}
		}

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[0] - F->tau2[0])));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau1[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau2[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau1[0])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau2[0])));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac * pow(prm.beta, 2);
	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac * pow(prm.beta, 2) * fac_ph;


	#if TEST_MODE
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif

	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){

		if( F->tau1[0] > F->tau2[0] ){
			F->flag=1;
		}
		else{
			F->flag=0;
		}
		#if TEST_MODE
		printf("  flag=%d\n", F->flag);
		#endif

		F->tau1[1] = F->tau1[0] + prm.beta;
		F->tau2[1] = F->tau2[0] + prm.beta;

// 		F->mat_Delta[0][0] = lambda;
		F->D.mat_M[0][0] = 1.0 / lambda;

		F->k++;
// 		K++;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void add_antisegment0(int sigma)
{
	Operators *F = &S[sigma];

	if(F->flag==1)  return;  // |0>

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}


	F->tau1[0] = rand_tau_beta(prm.beta);
	F->tau2[0] = rand_tau_beta(prm.beta);

	#if TEST_MODE
	printf("\nadd_antisegment0  (sigma=%d)\n", sigma);
	printf("  0  %8.5lf  %8.5lf\n", F->tau1[0], F->tau2[0]);
	#endif

	double lambda = G0_calc_interp(Delta[sigma], F->tau2[0], F->tau1[0]);


	double l_seg = F->tau2[0] - F->tau1[0];
	if(l_seg<0)  l_seg += prm.beta;

	double ln_b_fac = l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0
			double l_over = tau_overlap(F->tau1[0], F->tau2[0], S[s_ref[r]]);
			ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac += l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(F->tau1[0] - F->tau2[0]))
		                 - exp(-prm.w_ph * (prm.beta - fabs(F->tau1[0] - F->tau2[0])));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau1[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau1[0])))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau2[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau2[0])))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau1[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau1[0])))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau2[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau2[0])));
			}
		}

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[0] - F->tau2[0])));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau1[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau2[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau1[0])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau2[0])));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac * pow(prm.beta, 2);
	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac * pow(prm.beta, 2) * fac_ph;


	#if TEST_MODE
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif

	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){

		if( F->tau1[0] > F->tau2[0] ){
			F->flag=1;
		}
		else{
			F->flag=0;
		}
		#if TEST_MODE
		printf("  flag=%d\n", F->flag);
		#endif

		F->tau1[1] = F->tau1[0] + prm.beta;
		F->tau2[1] = F->tau2[0] + prm.beta;

// 		F->mat_Delta[0][0] = lambda;
		F->D.mat_M[0][0] = 1.0 / lambda;

		F->k++;
// 		K++;
		if(prob<0)  w_sign = -w_sign;
	}
}


static void add_tau(Operators &F, int i_tau1, double tau1, int i_tau2, double tau2)
{
	for(int i=F.k; i>i_tau1; i--)  F.tau1[i] = F.tau1[i-1];
	for(int i=F.k; i>i_tau2; i--)  F.tau2[i] = F.tau2[i-1];
	F.tau1[i_tau1] = tau1;
	F.tau2[i_tau2] = tau2;
	F.tau1[F.k+1] = F.tau1[0] + prm.beta;
	F.tau2[F.k+1] = F.tau2[0] + prm.beta;
}
static void rem_tau(Operators &F, int i_tau1, int i_tau2)
{
		for(int i=i_tau1; i<F.k-1; i++)  F.tau1[i] = F.tau1[i+1];
		for(int i=i_tau2; i<F.k-1; i++)  F.tau2[i] = F.tau2[i+1];
		F.tau1[F.k-1] = F.tau1[0] + prm.beta;
		F.tau2[F.k-1] = F.tau2[0] + prm.beta;
}

// *tau1 : annihilation
// *tau2 : creation
// tau_ins : creation
// i_tau_ins is stored
static int reject_create_seg(std::vector<double> &F_tau1, std::vector<double> &F_tau2, int k, int wind, double tau_ins, int &i_tau_ins, double &l_max)
{
// 	int wind = F->flag==1 ? 0 : 1;
	if( k ){
		// i_tau_ins = tau_order(F_tau2, k, tau_ins);
		i_tau_ins = tau_order(F_tau2, tau_ins);

		int i_tau1 = wind==0 ? i_tau_ins-1 : i_tau_ins;
		if( i_tau1 >= 0 && F_tau1[i_tau1] > tau_ins )  return 1;  // reject
		l_max = F_tau2[i_tau_ins] - tau_ins;
	}
	else{
		if( wind )  return 1;  // reject
		i_tau_ins = 0;
		l_max = prm.beta;
	}
	return 0;
}
static int reject_create_seg(Operators &F, double tau_ins, int &i_tau_ins, double &l_max)
{
	int wind = F.flag==1 ? 0 : 1;
	return reject_create_seg(F.tau1, F.tau2, F.k, wind, tau_ins, i_tau_ins, l_max);
}
/*
static int reject_create_seg(cond_op *F, double tau_ins, int &i_tau_ins, double &l_max)
{
	int wind = F->flag==1 ? 0 : 1;
	if( F->k ){
		i_tau_ins = tau_order(F->tau2, F->k, tau_ins);

		int i_tau1 = wind==0 ? i_tau_ins-1 : i_tau_ins;
		if( i_tau1 >= 0 && F->tau1[i_tau1] > tau_ins )  return 1;  // reject
		l_max = F->tau2[i_tau_ins] - tau_ins;
	}
	else{
		if( wind )  return 1;  // reject
		i_tau_ins = 0;
		l_max = prm.beta;
	}
	return 0;
}
*/

template <typename T>
void exchange_values(T &x, T &y)
{
	T temp = x;
	x = y;
	y = temp;
}

static void add_seg(int sigma, int anti)
{
	Operators *F = &S[sigma];

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	double tau1;  // for f-annihilation operator
	double tau2;  // for f-creation operator
	int i_tau1, i_tau2;
	double l_max, b_fac;
	int flag_change=0;

	{
		double tau_l, tau_r;
		int i_tau_l, i_tau_r;

		// double *F_tau_l = F->tau1, *F_tau_r = F->tau2;
		std::vector<double> &F_tau_l = F->tau1, &F_tau_r = F->tau2;

		int wind = F->flag==1 ? 0 : 1;
		if( anti ){
			exchange_values(F_tau_l, F_tau_r);
			wind ^= 1;
		}


		//
		// choose tau1 & tau2
		//
// 		tau2 = rand_tau_beta(prm.beta);  // f-creation
		tau_r = rand_tau_beta(prm.beta);  // f-creation

		#if TEST_MODE
		printf("\nadd_seg  (sigma=%d anti=%d)\n", sigma, anti);
		printf("  tau_r=%8.5lf\n", tau_r);
		#endif

// 		if( reject_create_seg(F, tau2, i_tau2, l_max) )  return;
// 		if( reject_create_seg(F->tau1, F->tau2, F->k, wind, tau2, i_tau2, l_max) )  return;
		if( reject_create_seg(F_tau_l, F_tau_r, F->k, wind, tau_r, i_tau_r, l_max) )  return;
		if( prm.UINF && anti==0 ){
			for(int r=0; r<N_S-1; r++){
				int i_tau_rtemp;
				double l_max_temp;
// 				if( reject_create_seg(&S[s_ref[r]], tau2, i_tau2_temp, l_max_temp) )  return;
				if( reject_create_seg(S[s_ref[r]], tau_r, i_tau_rtemp, l_max_temp) )  return;
				l_max = MIN(l_max, l_max_temp);
			}
		}

// 		i_tau1 = wind==0 ? i_tau2 : i_tau2 + 1;
		i_tau_l = wind==0 ? i_tau_r : i_tau_r + 1;
		double l_seg = rand_tau_l(l_max);
		tau_l = tau_r + l_seg;
		if( tau_l >= prm.beta ){
			tau_l -= prm.beta;
			flag_change = 1;
			i_tau_l = 0;
		}

		//
		// Boltzmann factor
		//
		double ln_b_fac = - l_seg * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
// 				double l_over = tau_overlap(tau2, tau1, S[s_ref[r]]);
				double l_over = tau_overlap(tau_r, tau_l, S[s_ref[r]]);
				ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
			}
		}
		if( anti )  ln_b_fac = -ln_b_fac;
		b_fac = exp(ln_b_fac);


		tau1 = tau_l;
		tau2 = tau_r;
		i_tau1 = i_tau_l;
		i_tau2 = i_tau_r;
		if( anti ){
			exchange_values(tau1, tau2);
			exchange_values(i_tau1, i_tau2);
		}
	}

	//
	// determinant
	//
	double Delta_tau1[F->k], Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}

	double diag = G0_calc_interp(Delta[sigma], tau2, tau1);
	// double lambda = add_lambda(*F, Delta_tau1, Delta_tau2, diag);
	double lambda = F->D.add_lambda(Delta_tau1, Delta_tau2, diag);


	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1);

	#if TEST_MODE
	printf("  tau1=%8.5lf  tau2=%8.5lf  l_max=%8.5lf\n", tau1, tau2, l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
		// add_tau(*F, i_tau1, tau1, i_tau2, tau2);
		add_tau(*F, i_tau1, tau1, i_tau2, tau2);

		F->flag = (F->flag + flag_change)%2;

		//
		// update mat_M
		//
		// add_mat_M(*F, i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);
		F->D.add_mat_M(i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);

		F->k++;
		F->tau1[F->k] = F->tau1[0] + prm.beta;
		F->tau2[F->k] = F->tau2[0] + prm.beta;

		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}
static void add_seg(int sigma)
{
	Operators *F = &S[sigma];

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	double tau1;  // for f-annihilation operator
	double tau2;  // for f-creation operator
	int i_tau1;
	int i_tau2;
	double l_max, l_seg;
	int flag_change=0;
	int wind = F->flag==1 ? 0 : 1;

	//
	// choose tau1 & tau2
	//
	tau2 = rand_tau_beta(prm.beta);  // f-creation
// 	i_tau2 = tau_order(F->tau2, F->k, tau2);

	#if TEST_MODE
	printf("\nadd_segment  (sigma=%d)\n", sigma);
	printf("  tau2=%8.5lf\n", tau2);
	#endif

	if( reject_create_seg(*F, tau2, i_tau2, l_max) )  return;
	if( prm.UINF ){
		for(int r=0; r<N_S-1; r++){
			int i_tau2_temp;
			double l_max_temp;
			if( reject_create_seg(S[s_ref[r]], tau2, i_tau2_temp, l_max_temp) )  return;
			l_max = MIN(l_max, l_max_temp);
		}
	}

	i_tau1 = wind==0 ? i_tau2 : i_tau2 + 1;
	l_seg = rand_tau_l(l_max);
	tau1 = tau2 + l_seg;
	if( tau1 >= prm.beta ){
		tau1 -= prm.beta;
		flag_change = 1;
		i_tau1 = 0;
	}

	//
	// Boltzmann factor
	//
	double ln_b_fac = - l_seg * prm.ef[sigma];
	if( prm.UINF==0 ){
		for(int r=0; r<N_S-1; r++){
			double l_over = tau_overlap(tau2, tau1, S[s_ref[r]]);
			ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
		}
	}
	double b_fac = exp(ln_b_fac);

	//
	// determinant
	//
	double Delta_tau1[F->k], Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}

	double diag = G0_calc_interp(Delta[sigma], tau2, tau1);
	// double lambda = add_lambda(*F, Delta_tau1, Delta_tau2, diag);
	double lambda = F->D.add_lambda(Delta_tau1, Delta_tau2, diag);


	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1);

	#if TEST_MODE
	printf("  tau1=%8.5lf  l_max=%8.5lf\n", tau1, l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
		add_tau(*F, i_tau1, tau1, i_tau2, tau2);

		F->flag = (F->flag + flag_change)%2;

		//
		// update mat_M
		//
		// add_mat_M(*F, i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);
		F->D.add_mat_M(i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);

		F->k++;
		F->tau1[F->k] = F->tau1[0] + prm.beta;
		F->tau2[F->k] = F->tau2[0] + prm.beta;

		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}

static void rem_seg(int sigma, int anti)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		int sigma2 = rand_int(N_S);
		if( sigma == sigma2 )  state0_change(sigma);
		else  state0_change_2(sigma, sigma2);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nrem_seg  (sigma=%d anti=%d)\n", sigma, anti);
	#endif

	int flag_change=0;
	int i_tau1, i_tau2;
	double l_max, b_fac;
	{
		// double *F_tau_l = F->tau1, *F_tau_r = F->tau2;
		std::vector<double> *F_tau_l = &(F->tau1), *F_tau_r = &(F->tau2);
		int wind = F->flag==1 ? 0 : 1;
		if( anti ){
			exchange_values(F_tau_l, F_tau_r);
			wind ^= 1;
		}

		//
		// choose tau1 & tau2
		//
		int i_tau_r = rand_int(F->k);
		int i_tau_l = wind==0 ? i_tau_r : i_tau_r + 1;  // may be i_tau_l==F->k

		double l_seg = (*F_tau_l)[i_tau_l] - (*F_tau_r)[i_tau_r];
		if( i_tau_l == F->k ){
			i_tau_l = 0;
			flag_change = 1;
		}

		l_max = prm.beta;
// 		if( F->k > 1 )  l_max = F->tau2[i_tau2+1] - F->tau2[i_tau2];
		if( F->k > 1 )  l_max = (*F_tau_r)[i_tau_r+1] - (*F_tau_r)[i_tau_r];

		if( prm.UINF && anti==0 ){
			for(int r=0; r<N_S-1; r++){
				int i_tau_rtemp;
				double l_max_temp;
// 				double tau2 = F->tau2[i_tau2];
				double tau_r = (*F_tau_r)[i_tau_r];
				reject_create_seg(S[s_ref[r]], tau_r, i_tau_rtemp, l_max_temp);
				l_max = MIN(l_max, l_max_temp);
			}
		}
		if( prm.UINF && anti ){
// 			printf("     tau_r=%lf\n", (*F_tau_r)[i_tau_r]);
// 			printf("     tau_l=%lf\n", (*F_tau_l)[i_tau_l]);
			for(int r=0; r<N_S-1; r++){
				double l_over = tau_overlap((*F_tau_r)[i_tau_r], (*F_tau_l)[i_tau_l], S[s_ref[r]]);
// 				printf("     sigma=%d  l_over=%lf\n", s_ref[r], l_over);
				if( l_over )  return;
			}
		}


		//
		// Boltzmann factor
		//
		double ln_b_fac = l_seg * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
// 				double l_over = tau_overlap(F->tau2[i_tau2], F->tau1[i_tau1], S[s_ref[r]]);
				double l_over = tau_overlap((*F_tau_r)[i_tau_r], (*F_tau_l)[i_tau_l], S[s_ref[r]]);
				ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
			}
		}
		if( anti )  ln_b_fac = -ln_b_fac;
		b_fac = exp(ln_b_fac);


		i_tau1 = i_tau_l;
		i_tau2 = i_tau_r;
		if( anti ){
			exchange_values(i_tau1, i_tau2);
		}
	}

	double lambda = F->D.mat_M[i_tau1][i_tau2];

	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max);

	#if TEST_MODE
	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void remove_tau(Operators &F, int i_tau1, int i_tau2);
		rem_tau(*F, i_tau1, i_tau2);

		F->flag = (F->flag + flag_change)%2;

		//
		// update mat_M : delete row [i_tau1] & column [i_tau2]
		//
		// remove_mat_M(*F, i_tau1, i_tau2, lambda);
		F->D.remove_mat_M(i_tau1, i_tau2, lambda);

		F->k--;
		F->tau1[F->k] = F->tau1[0] + prm.beta;
		F->tau2[F->k] = F->tau2[0] + prm.beta;

		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}
static void rem_seg(int sigma)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		int sigma2 = rand_int(N_S);
		if( sigma == sigma2 )  state0_change(sigma);
		else  state0_change_2(sigma, sigma2);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nremove_segment  (sigma=%d)\n", sigma);
	#endif

	//
	// choose tau1 & tau2
	//
	int flag_change=0;
	int wind = F->flag==1 ? 0 : 1;

	int i_tau2 = rand_int(F->k);
	int i_tau1 = wind==0 ? i_tau2 : i_tau2 + 1;

	double l_max = prm.beta;
	if( F->k > 1 )  l_max = F->tau2[i_tau2+1] - F->tau2[i_tau2];

	if( prm.UINF ){
		for(int r=0; r<N_S-1; r++){
			int i_tau2_temp;
			double l_max_temp;
			double tau2 = F->tau2[i_tau2];
			// reject_create_seg(&S[s_ref[r]], tau2, i_tau2_temp, l_max_temp);
			reject_create_seg(S[s_ref[r]], tau2, i_tau2_temp, l_max_temp);
			l_max = MIN(l_max, l_max_temp);
		}
	}

	//
	// Boltzmann factor
	//
	double l_seg = F->tau1[i_tau1] - F->tau2[i_tau2];

	if( i_tau1 == F->k ){
		i_tau1 = 0;
		flag_change = 1;
	}

	double ln_b_fac = l_seg * prm.ef[sigma];
	if( prm.UINF==0 ){
		for(int r=0; r<N_S-1; r++){
			double l_over = tau_overlap(F->tau2[i_tau2], F->tau1[i_tau1], S[s_ref[r]]);
			ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
		}
	}
	double b_fac = exp(ln_b_fac);


	double lambda = F->D.mat_M[i_tau1][i_tau2];

	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max);

	#if TEST_MODE
	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void remove_tau(Operators &F, int i_tau1, int i_tau2);
		rem_tau(*F, i_tau1, i_tau2);

		F->flag = (F->flag + flag_change)%2;

		//
		// update mat_M : delete row [i_tau1] & column [i_tau2]
		//
		// remove_mat_M(*F, i_tau1, i_tau2, lambda);
		F->D.remove_mat_M(i_tau1, i_tau2, lambda);

		F->k--;
		F->tau1[F->k] = F->tau1[0] + prm.beta;
		F->tau2[F->k] = F->tau2[0] + prm.beta;

		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}

static void add_segment(int sigma)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		add_segment0(sigma);
// 		add_0(sigma);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	double tau1;  // for f-annihilation operator
	double tau2;  // for f-creation operator
	int i_tau1;
	int i_tau2;
	double l_max, l_seg;
	int flag_change=0;

	//
	// choose tau1 & tau2
	//

	tau1 = rand_tau_beta(prm.beta);  // f-annihilation
	// i_tau1 = tau_order(F->tau1, F->k, tau1);
	i_tau1 = tau_order(F->tau1, tau1);

	#if TEST_MODE
	printf("\nadd_segment  (sigma=%d)\n", sigma);
	printf("  i_tau1=%d  tau1=%8.5lf\n", i_tau1, tau1);
	#endif

	if(F->flag){  // flag==1: tau1[i] > tau2[i]
		if( tau1 > F->tau2[i_tau1] ){  // cut a segment
			return;
		}
		else{  // add a segment
			if( i_tau1==0 ){
				l_max = tau1 + prm.beta - F->tau1[F->k-1];
				l_seg = rand_tau_l(l_max);
				tau2 = tau1 - l_seg;
				if(tau2<0){
					tau2 += prm.beta;
					i_tau2 = F->k;
					flag_change = 1;
				}
				else{
					i_tau2 = i_tau1;
				}
			}
			else{
				l_max = tau1 - F->tau1[i_tau1-1];
				l_seg = rand_tau_l(l_max);
				tau2 = tau1 - l_seg;
				i_tau2 = i_tau1;
			}
		}
	}
	else{  // flag==0: tau1[i] < tau2[i]
		if( i_tau1==0 ||  tau1 > F->tau2[i_tau1-1] ){  // cut a segment
			return;
		}
		else{  // add a segment
			i_tau2 = i_tau1 - 1;
			l_max = tau1 - F->tau1[i_tau1-1];
			l_seg = rand_tau_l(l_max);
			tau2 = tau1 - l_seg;
		}
	}

	//
	// Boltzmann factor
	//
	double ln_b_fac = - l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0

			double l_over = tau_overlap(tau2, tau1, S[s_ref[r]]);
			ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac -= l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	//
	// determinant
	//
	double Delta_tau1[F->k], Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}

	double diag = G0_calc_interp(Delta[sigma], tau2, tau1);
	// double lambda = add_lambda(*F, Delta_tau1, Delta_tau2, diag);
	double lambda = F->D.add_lambda(Delta_tau1, Delta_tau2, diag);


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(tau1-tau2))
		                 - exp(-prm.w_ph * (prm.beta - fabs(tau1-tau2)));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - tau1))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau1)))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - tau2))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau2)))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - tau1))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau1)))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - tau2))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau2)));
			}
		}

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(tau1-tau2)));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau1)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau2)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau1)))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau2)));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1);
	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1) * fac_ph;

	#if TEST_MODE
 	printf("  i_tau2=%d  tau2=%8.5lf  l_max=%8.5lf\n", i_tau2, tau2, l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void add_tau(Operators &F, int, double, int, double);
		add_tau(*F, i_tau1, tau1, i_tau2, tau2);

		F->flag = (F->flag + flag_change)%2;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<=F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif

		//
		// update mat_M
		//
		// add_mat_M(*F, i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);
		F->D.add_mat_M(i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);


		F->k++;
// 		K++;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void add_antisegment(int sigma)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		add_antisegment0(sigma);
// 		add_0(sigma);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	double tau1;  // for f-annihilation operator
	double tau2;  // for f-creation operator
	int i_tau1;
	int i_tau2;
	double l_max, l_seg;
	int flag_change=0;

	//
	// choose tau1 & tau2
	//
	tau1 = rand_tau_beta(prm.beta);  // f-annihilation
	// i_tau1 = tau_order(F->tau1, F->k, tau1);
	i_tau1 = tau_order(F->tau1, tau1);

	#if TEST_MODE
	printf("\nadd_antisegment  (sigma=%d)\n", sigma);
	printf("  i_tau1=%d  tau1=%8.5lf\n", i_tau1, tau1);
	#endif

	if(F->flag){  // flag==1: tau1[i] > tau2[i]
		if( tau1 > F->tau2[i_tau1] ){  // cut a segment
			i_tau2 = i_tau1 + 1;
			l_max = F->tau1[i_tau1] - tau1;
			l_seg = rand_tau_l(l_max);
			tau2 = tau1 + l_seg;
		}
		else{  // add a segment
			return;
		}
	}
	else{  // flag==0: tau1[i] < tau2[i]
		if( i_tau1==0 ){  // cut a segment
			l_max = F->tau1[i_tau1] - tau1;
			l_seg = rand_tau_l(l_max);
			tau2 = tau1 + l_seg;
			i_tau2 = i_tau1;
		}
		else if( tau1 > F->tau2[i_tau1-1] ){  // cut a segment
			l_max = F->tau1[i_tau1] - tau1;
			l_seg = rand_tau_l(l_max);
			tau2 = tau1 + l_seg;
			if( tau2>=prm.beta ){
				tau2 -= prm.beta;
				i_tau2 = 0;
				flag_change = 1;
			}
			else{
				i_tau2 = i_tau1;
			}
		}
		else{  // add a segment
			return;
		}
	}

	//
	// Boltzmann factor
	//
	double ln_b_fac = l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0

			double l_over = tau_overlap(tau1, tau2, S[s_ref[r]]);
			ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac += l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	//
	// determinant
	//
	double Delta_tau1[F->k], Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}

	double diag = G0_calc_interp(Delta[sigma], tau2, tau1);
	// double lambda = add_lambda(*F, Delta_tau1, Delta_tau2, diag);
	double lambda = F->D.add_lambda(Delta_tau1, Delta_tau2, diag);


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(tau1-tau2))
		                 - exp(-prm.w_ph * (prm.beta - fabs(tau1-tau2)));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - tau1))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau1)))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - tau2))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau2)))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - tau1))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau1)))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - tau2))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau2)));
			}
		}

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(tau1-tau2)));

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau1)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau2)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau1)))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau2)));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1);
	double prob = sgn(tau1-tau2) * lambda * b_fac * prm.beta * l_max / (double)(F->k+1) * fac_ph;

	#if TEST_MODE
	printf("  i_tau2=%d  tau2=%8.5lf  l_max=%8.5lf\n", i_tau2, tau2, l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void add_tau(Operators &F, int, double, int, double);
		add_tau(*F, i_tau1, tau1, i_tau2, tau2);

		F->flag = (F->flag + flag_change)%2;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<=F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif

		//
		// update mat_M
		//
		// add_mat_M(*F, i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);
		F->D.add_mat_M(i_tau1, i_tau2, Delta_tau1, Delta_tau2, lambda);


		F->k++;
// 		K++;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void remove_segment0(int sigma)
{
	Operators *F = &S[sigma];

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nremove_segment0  (sigma=%d)\n", sigma);
	#endif


	double l_seg = F->tau1[0] - F->tau2[0];
	if(l_seg<0)  l_seg += prm.beta;

	double ln_b_fac = l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0
			double l_over = tau_overlap(F->tau2[0], F->tau1[0], S[s_ref[r]]);
			ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac += l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);

	double lambda = F->D.mat_M[0][0];


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(F->tau1[0] - F->tau2[0]))
		                 - exp(-prm.w_ph * (prm.beta - fabs(F->tau1[0] - F->tau2[0])));

		for(int r=0; r<N_S-1; r++){
			for(int i=0; i<S[s_ref[r]].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s_ref[r]].tau1[i] - F->tau1[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau1[i] - F->tau1[0])))
				           - exp(-prm.w_ph * fabs(S[s_ref[r]].tau1[i] - F->tau2[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau1[i] - F->tau2[0])))
				           - exp(-prm.w_ph * fabs(S[s_ref[r]].tau2[i] - F->tau1[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau2[i] - F->tau1[0])))
				           + exp(-prm.w_ph * fabs(S[s_ref[r]].tau2[i] - F->tau2[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau2[i] - F->tau2[0])));
			}
		}

		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[0] - F->tau2[0])));

		for(int r=0; r<N_S-1; r++){
			for(int i=0; i<S[s_ref[r]].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau1[i] - F->tau1[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau1[i] - F->tau2[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau2[i] - F->tau1[0])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau2[i] - F->tau2[0])));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac / pow(prm.beta, 2);
	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac / pow(prm.beta, 2) * fac_ph;

	#if TEST_MODE
// 	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){

		F->flag=1;  // |0>

		F->k=0;
// 		K--;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void remove_antisegment0(int sigma)
{
	Operators *F = &S[sigma];

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nremove_antisegment0  (sigma=%d)\n", sigma);
	#endif


	double l_seg = F->tau2[0] - F->tau1[0];
	if(l_seg<0)  l_seg += prm.beta;

	double ln_b_fac = -l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0
			double l_over = tau_overlap(F->tau1[0], F->tau2[0], S[s_ref[r]]);
			ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac -= l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);

	double lambda = F->D.mat_M[0][0];


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph)
		                 - exp(-prm.w_ph * fabs(F->tau1[0] - F->tau2[0]))
		                 - exp(-prm.w_ph * (prm.beta - fabs(F->tau1[0] - F->tau2[0])));

		for(int r=0; r<N_S-1; r++){
			for(int i=0; i<S[s_ref[r]].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s_ref[r]].tau1[i] - F->tau1[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau1[i] - F->tau1[0])))
				           - exp(-prm.w_ph * fabs(S[s_ref[r]].tau1[i] - F->tau2[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau1[i] - F->tau2[0])))
				           - exp(-prm.w_ph * fabs(S[s_ref[r]].tau2[i] - F->tau1[0]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau2[i] - F->tau1[0])))
				           + exp(-prm.w_ph * fabs(S[s_ref[r]].tau2[i] - F->tau2[0]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s_ref[r]].tau2[i] - F->tau2[0])));
			}
		}

		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5)
		                 - cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[0] - F->tau2[0])));

		for(int r=0; r<N_S-1; r++){
			for(int i=0; i<S[s_ref[r]].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau1[i] - F->tau1[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau1[i] - F->tau2[0])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau2[i] - F->tau1[0])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s_ref[r]].tau2[i] - F->tau2[0])));
			}
		}

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac / pow(prm.beta, 2);
	double prob = sgn(F->tau1[0]-F->tau2[0]) * lambda * b_fac / pow(prm.beta, 2) * fac_ph;

	#if TEST_MODE
// 	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){

		F->flag=0;  // |1>

		F->k=0;
// 		K--;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void remove_tau(Operators &F, int i_tau1, int i_tau2)
{
		for(int i=i_tau1; i<F.k-1; i++)  F.tau1[i] = F.tau1[i+1];
		for(int i=i_tau2; i<F.k-1; i++)  F.tau2[i] = F.tau2[i+1];
		F.tau1[F.k-1] = F.tau1[0] + prm.beta;
		F.tau2[F.k-1] = F.tau2[0] + prm.beta;
}


static void remove_segment(int sigma)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		int sigma2 = rand_int(N_S);
		if( sigma == sigma2 )  state0_change(sigma);
		else  state0_change_2(sigma, sigma2);
		return;
	}
	else if(F->k==1){
		remove_segment0(sigma);
// 		remove_0(sigma);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nremove_segment  (sigma=%d)\n", sigma);
	#endif

	//
	// choose tau1 & tau2
	//
	int i_tau1 = rand_int(F->k);
	int i_tau2;
	int flag_change=0;

	if( F->flag ){  // flag==1: tau1[i] > tau2[i]
		i_tau2 = i_tau1;
	}
	else{
		if( i_tau1==0 ){
			i_tau2 = F->k-1;
			flag_change = 1;
		}
		else  i_tau2 = i_tau1 - 1;
	}


	//
	// Boltzmann factor
	//
	double l_seg = F->tau1[i_tau1] - F->tau2[i_tau2];
	if(l_seg < 0)  l_seg += prm.beta;

	double ln_b_fac = l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0

			double l_over = tau_overlap(F->tau2[i_tau2], F->tau1[i_tau1], S[s_ref[r]]);
			ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac += l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	double l_max;
	if(i_tau1==0){
		l_max = prm.beta + F->tau1[i_tau1] - F->tau1[F->k-1];
	}
	else{
		l_max = F->tau1[i_tau1] - F->tau1[i_tau1-1];
	}

	double lambda = F->D.mat_M[i_tau1][i_tau2];


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph);

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau1[i_tau1]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau2[i_tau2]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau2[i_tau2])))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau1[i_tau1]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau1[i_tau1])))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau2[i_tau2]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= 2.0 * (1.0 + exp(-prm.beta*prm.w_ph));
		ln_fac_ph += exp(-prm.w_ph * fabs(F->tau1[i_tau1] - F->tau2[i_tau2]))
		           + exp(-prm.w_ph * (prm.beta - fabs(F->tau1[i_tau1] - F->tau2[i_tau2])));

		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5);

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau2[i_tau2])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau1[i_tau1])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= 2.0 * cosh(prm.w_ph * prm.beta*0.5);
		ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[i_tau1] - F->tau2[i_tau2])));

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max);
	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max) * fac_ph;

	#if TEST_MODE
	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void remove_tau(Operators &F, int i_tau1, int i_tau2);
		remove_tau(*F, i_tau1, i_tau2);

		F->flag = (F->flag + flag_change)%2;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif

		//
		// update mat_M : delete row [i_tau1] & column [i_tau2]
		//
		// remove_mat_M(*F, i_tau1, i_tau2, lambda);
		F->D.remove_mat_M(i_tau1, i_tau2, lambda);


		F->k--;
// 		K--;

// 		if(F->k==0)  F->flag=1;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void remove_antisegment(int sigma)
{
	Operators *F = &S[sigma];

	if(F->k==0){
		int sigma2 = rand_int(N_S);
		if( sigma == sigma2 )  state0_change(sigma);
		else  state0_change_2(sigma, sigma2);
		return;
	}
	else if(F->k==1){
		remove_antisegment0(sigma);
// 		remove_0(sigma);
		return;
	}

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if TEST_MODE
	printf("\nremove_antisegment  (sigma=%d)\n", sigma);
	#endif

	//
	// choose tau1 & tau2
	//
	int i_tau1 = rand_int(F->k);
	int i_tau2;
	int flag_change=0;

	if( F->flag ){  // flag==1: tau1[i] > tau2[i]
		if( i_tau1==F->k-1 ){
			i_tau2 = 0;
			flag_change = 1;
		}
		else  i_tau2 = i_tau1 + 1;
	}
	else{
		i_tau2 = i_tau1;
	}


	//
	// Boltzmann factor
	//
	double l_seg = F->tau2[i_tau2] - F->tau1[i_tau1];
	if(l_seg < 0)  l_seg += prm.beta;

	double ln_b_fac = -l_seg * prm.ef[sigma];

	for(int r=0; r<N_S-1; r++){
		if(S[s_ref[r]].k){  // !=0

			double l_over = tau_overlap(F->tau1[i_tau1], F->tau2[i_tau2], S[s_ref[r]]);
			ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
		}
		else{
			if(S[s_ref[r]].flag==0){  // |1>
				ln_b_fac -= l_seg * prm.U[sigma][s_ref[r]];
			}
		}
	}

	double b_fac = exp(ln_b_fac);


	double l_max = F->tau1[i_tau1+1] - F->tau1[i_tau1];

	double lambda = F->D.mat_M[i_tau1][i_tau2];


	//
	// phonon part
	//
	double fac_ph=1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 1.0 + exp(-prm.beta * prm.w_ph);

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau1[i_tau1]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau2[i_tau2]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau2[i_tau2])))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau1[i_tau1]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau1[i_tau1])))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau2[i_tau2]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= 2.0 * (1.0 + exp(-prm.beta*prm.w_ph));
		ln_fac_ph += exp(-prm.w_ph * fabs(F->tau1[i_tau1] - F->tau2[i_tau2]))
		           + exp(-prm.w_ph * (prm.beta - fabs(F->tau1[i_tau1] - F->tau2[i_tau2])));

		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = cosh(prm.w_ph * prm.beta*0.5);

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau2[i_tau2])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau1[i_tau1])))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= 2.0 * cosh(prm.w_ph * prm.beta*0.5);
		ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[i_tau1] - F->tau2[i_tau2])));

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max);
	double prob = sgn(F->tau1[i_tau1]-F->tau2[i_tau2]) * lambda * b_fac * (double)(F->k) / (prm.beta * l_max) * fac_ph;

	#if TEST_MODE
	printf(" %d %8.5lf   %d %8.5lf\n", i_tau1, F->tau1[i_tau1], i_tau2, F->tau2[i_tau2]);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_seg, ACCPT_sub.n_reject_seg) ){
		//
		// update tau-array
		//
// 		void remove_tau(Operators &F, int i_tau1, int i_tau2);
		remove_tau(*F, i_tau1, i_tau2);

		F->flag = (F->flag + flag_change)%2;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif

		//
		// update mat_M : delete row [i_tau1] & column [i_tau2]
		//
		// remove_mat_M(*F, i_tau1, i_tau2, lambda);
		F->D.remove_mat_M(i_tau1, i_tau2, lambda);


		F->k--;
// 		K--;

// 		if(F->k==0)  F->flag=0;
		if(prob<0)  w_sign = -w_sign;
	}
}

//
// transition between |0> and |1>
//

static void state0_change(int sigma)
{
	if( S[sigma].k )  return;

	int s_ref[N_S];  // only N_S-1 elements are used
// 	for(int r=0; r<N_S; r++)  s_ref[r] = r;
// 	for(int r=sigma; r<N_S; r++)  s_ref[r] ++;
// 	for(int r=0; r<N_S-1; r++){
// 		if(r<sigma)  s_ref[r] = r;
// 		else         s_ref[r] = r+1;
// 	}
	for(int r=0; r<N_S; r++)  s_ref[r] = r;
	s_ref[sigma] = N_S-1;

	double ln_b_fac = prm.beta * prm.ef[sigma];
	if( prm.UINF ){
		for(int r=0; r<N_S-1; r++){
			if( S[s_ref[r]].k || S[s_ref[r]].flag==0 )  return;
		}
	}
	else{
		for(int r=0; r<N_S-1; r++){
			ln_b_fac += tau_length(S[s_ref[r]]) * prm.U[sigma][s_ref[r]];
		}
	}

	if( S[sigma].flag )  ln_b_fac = -ln_b_fac;  // |0> -> |1>
	double prob = exp(ln_b_fac);

	unsigned long n_accept_temp=0, n_reject_temp=0;

// 	if( metropolis(prob, n_accept_temp, n_reject_temp) ){
	if( metropolis_abs(fabs(prob), n_accept_temp, n_reject_temp) ){

		S[sigma].flag = (S[sigma].flag + 1)%2;
		if(prob<0)  w_sign = -w_sign;
	}
}

static void state0_change_2(int sigma1, int sigma2)
{
	if( S[sigma1].k || S[sigma2].k || sigma1==sigma2 )  return;

	int s_ref[N_S];  // only N_S-2 elements are used
	if( sigma1 > sigma2 )  std::swap(sigma1, sigma2);
	for(int r=0; r<N_S-2; r++){
		s_ref[r] = r;
		if(s_ref[r] >= sigma1 )  s_ref[r]++;
		if(s_ref[r] >= sigma2 )  s_ref[r]++;
	}

	// check
// 	{
// 		printf("%d %d  ", sigma1, sigma2);
// 		for(int r=0; r<N_S-2; r++){
// 			printf(" %d", s_ref[r]);
// 		}
// 		printf("\n");
//
// 		int flag;
// 		for(int r=0; r<N_S; r++)  flag ^= 1 << r;
// 		for(int r=0; r<N_S-2; r++)  flag ^= 1 << s_ref[r];
// 		flag ^= 1 << sigma1;
// 		flag ^= 1 << sigma2;
// 		if( flag ){
// 			printf(" ERROR: stat0_change_2\n");
// 			exit(0);
// 		}
// 	}

	double ln_b_fac1 = prm.beta * prm.ef[sigma1];
	double ln_b_fac2 = prm.beta * prm.ef[sigma2];
	if( prm.UINF ){
		if( S[sigma1].flag == S[sigma2].flag )  return;
		for(int r=0; r<N_S-2; r++){
			if( S[s_ref[r]].k || S[s_ref[r]].flag==0 )  return;
		}
	}
	else{
		for(int r=0; r<N_S-2; r++){
			double l = tau_length(S[s_ref[r]]);
			ln_b_fac1 += tau_length(S[s_ref[r]]) * prm.U[sigma1][s_ref[r]];
			ln_b_fac2 += tau_length(S[s_ref[r]]) * prm.U[sigma2][s_ref[r]];
		}
	}

	if( S[sigma1].flag )  ln_b_fac1 = -ln_b_fac1;  // |0> -> |1>
	if( S[sigma2].flag )  ln_b_fac2 = -ln_b_fac2;

	double ln_u = 0;
	if( S[sigma1].flag == S[sigma2].flag )  ln_u = prm.beta * prm.U[sigma1][sigma2];
	if( S[sigma1].flag )  ln_u = -ln_u;  // |0> -> |1>

	double prob = exp(ln_b_fac1 + ln_b_fac2 + ln_u);

	unsigned long n_accept_temp=0, n_reject_temp=0;

	if( metropolis_abs(fabs(prob), n_accept_temp, n_reject_temp) ){

		S[sigma1].flag = (S[sigma1].flag + 1)%2;
		S[sigma2].flag = (S[sigma2].flag + 1)%2;
		if(prob<0)  w_sign = -w_sign;
	}
}

/*
//
// transition between |0> and |1>
//
static void state0_change(int sigma)
{
// 	if(S[sigma].k)  return;

	int s_ref[N_S];  // only N_S-1 elements are used
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

	#if N_S==2
	{
		if(S[s_ref[0]].k==0){
			if(S[s_ref[0]].flag == S[sigma].flag){
				if(rand_int(2)){
					state0_change_empty_double();
					return;
				}
			}
			else{
				if(rand_int(2)){
					state0_exchange(sigma, s_ref[0]);
					return;
				}
			}
		}
	}
	#else // N_S==2
	{
		// list of components with which 'flag' can be exchanged
		int n_exc=0, s_exc[N_S];  // only N_S-1 elements are used
		for(int r=0; r<N_S-1; r++){
			if(S[s_ref[r]].k==0){
				if(S[s_ref[r]].flag != S[sigma].flag)   s_exc[n_exc++] = s_ref[r];
			}
		}

		int i = rand_int(n_exc+1);
		if(i < n_exc){
			state0_exchange(sigma, s_exc[i]);
			return;
		}
	}
	#endif // N_S==2


	double ln_b_fac = prm.beta * prm.ef[sigma];
	for(int r=0; r<N_S-1; r++){
		ln_b_fac += tau_length(S[s_ref[r]]) * prm.U[sigma][s_ref[r]];
	}

	double prob;
	if(S[sigma].flag)  prob = exp(-ln_b_fac);  // |0> -> |1>
	else               prob = exp( ln_b_fac);  // |1> -> |0>


	unsigned long n_accept_temp=0, n_reject_temp=0;

// 	if( metropolis(prob, n_accept_temp, n_reject_temp) ){
	if( metropolis_abs(fabs(prob), n_accept_temp, n_reject_temp) ){

		S[sigma].flag = (S[sigma].flag + 1)%2;
		if(prob<0)  w_sign = -w_sign;
	}
}
*/
// this function is called only when N_S >= 2
static void state0_exchange(int sigma1, int sigma2)
{

	int s_ref[N_S];  // only N_S-2 elements are used
	for(int r=0; r<N_S-2; r++){
		if(r<sigma1){
			if(r<sigma2)  s_ref[r] = r;
			else          s_ref[r] = r+1;
		}
		else{
			if(r<sigma2)  s_ref[r] = r+1;
			else          s_ref[r] = r+2;
		}
	}

	double ln_b_fac = prm.beta * (prm.ef[sigma1] - prm.ef[sigma2]);
	for(int r=0; r<N_S-2; r++){
		ln_b_fac += tau_length(S[s_ref[r]]) * (prm.U[sigma1][s_ref[r]] - prm.U[sigma2][s_ref[r]]);
	}

	double prob;
	if(S[sigma1].flag)  prob = exp(-ln_b_fac);  // sigma1: |0> -> |1>,  sigma2: |1> -> |0>
	else                prob = exp( ln_b_fac);  // sigma1: |1> -> |0>,  sigma2: |0> -> |1>


	unsigned long n_accept_temp=0, n_reject_temp=0;

// 	if( metropolis(prob, n_accept_temp, n_reject_temp) ){
	if( metropolis_abs(fabs(prob), n_accept_temp, n_reject_temp) ){

		S[sigma1].flag = (S[sigma1].flag + 1)%2;
		S[sigma2].flag = (S[sigma2].flag + 1)%2;
		if(prob<0)  w_sign = -w_sign;
	}
}

// this function is called only when N_S==2
static void state0_change_empty_double()
{
	double ln_fac_b = prm.beta * (prm.ef[0] + prm.ef[1] + prm.U[0][1]);

	double prob;
	if(S[0].flag)  prob = exp(-ln_fac_b);  // empty -> double
	else           prob = exp( ln_fac_b);  // double -> empty


	unsigned long n_accept_temp=0, n_reject_temp=0;

// 	if( metropolis(prob, n_accept_temp, n_reject_temp) ){
	if( metropolis_abs(fabs(prob), n_accept_temp, n_reject_temp) ){

		S[0].flag = (S[0].flag + 1)%2;
		S[1].flag = (S[1].flag + 1)%2;
		if(prob<0)  w_sign = -w_sign;
	}
}


static void shift_tau1_rand(int sigma)
{
	if(S[sigma].k)  shift_tau1(sigma, rand_int(S[sigma].k));;
}

static void shift_tau1_all(int sigma)
{
	if(S[sigma].flag){  // flag==1: tau1[i] > tau2[i],  execute in ascending order
		for(int i=0; i<S[sigma].k; i++){
			shift_tau1(sigma, i);
		}
	}
	else{  // flag==0: tau1[i] < tau2[i],  execute in descending order
		for(int i=S[sigma].k-1; i>=0; i--){
			shift_tau1(sigma, i);
		}
	}

	#if TEST_MODE
	printf("\n");
	for(int i=0; i<S[sigma].k; i++){
		printf("%3d  %8.5lf  %8.5lf\n", i, S[sigma].tau1[i], S[sigma].tau2[i]);
	}
	printf("  flag=%d\n", S[sigma].flag);
	#endif
}

//
// shift tau1 (f-annihilation op., c-creation op.)
//
static void shift_tau1(int sigma, int i_tau1)
{
	Operators *F = &S[sigma];

	int s_ref[N_S-1];
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

// 	if(F->k==0)  return;
//
// 	int i_tau1 = rand_int(F->k);

// 	double l_over, l_total;
// 	double b_fac;  // Boltzmann factor
	double l_max, l_dif, tau1;
	int flag_change=0;
// 	int flag_direction;  // 0: segment shrinks,  1: segment extends

	if( F->flag ){  // flag==1: tau1[i] > tau2[i]
		l_max = F->tau2[i_tau1+1] - F->tau2[i_tau1];
		if( prm.UINF ){
			for(int r=0; r<N_S-1; r++){
				// cond_op *p_S = &S[s_ref[r]];
				Operators *p_S = &S[s_ref[r]];
				if( p_S->k ){
					// int i_tau2_temp = tau_order(p_S->tau2, p_S->k, F->tau1[i_tau1]);
					int i_tau2_temp = tau_order(p_S->tau2, F->tau1[i_tau1]);
					double l_max_temp = p_S->tau2[i_tau2_temp] - F->tau2[i_tau1];
					l_max = MIN(l_max, l_max_temp);
				}
			}
		}

		tau1 = F->tau2[i_tau1] + rand_tau_l(l_max);
		l_dif = tau1 - F->tau1[i_tau1];

		if( tau1 > prm.beta ){
			flag_change = 1;
			tau1 -= prm.beta;
		}
	}
	else{  // flag==0: tau1[i] < tau2[i]
		int i_tau1m = (i_tau1 - 1 + F->k) % F->k;
		l_max = F->tau2[i_tau1] - F->tau2[i_tau1m];
		if(i_tau1==0)  l_max += prm.beta;

		double tau2_temp = F->tau2[i_tau1];
		if( prm.UINF ){
			for(int r=0; r<N_S-1; r++){
				// cond_op *p_S = &S[s_ref[r]];
				Operators *p_S = &S[s_ref[r]];
				if( p_S->k ){
					// int i_tau2_temp = tau_order(p_S->tau2, p_S->k, F->tau1[i_tau1]);
					int i_tau2_temp = tau_order(p_S->tau2, F->tau1[i_tau1]);
					double l_max_temp = p_S->tau2[i_tau2_temp] - F->tau2[i_tau1m];
					if( l_max_temp < 0 )  l_max_temp += prm.beta;
					if( l_max > l_max_temp ){
						l_max = l_max_temp;
						tau2_temp = p_S->tau2[i_tau2_temp%p_S->k];
					}
				}
			}
		}

// 		tau1 = F->tau2[i_tau1] - rand_tau_l(l_max);
		tau1 = tau2_temp - rand_tau_l(l_max);
		l_dif = tau1 - F->tau1[i_tau1];

		if( tau1 < 0 ){
			flag_change = 1;
			tau1 += prm.beta;
		}
	}

	double ln_b_fac;
	if( l_dif > 0){ // segment extends
		ln_b_fac = -l_dif * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
				if(S[s_ref[r]].k){
					double l_over = tau_overlap(F->tau1[i_tau1], tau1, S[s_ref[r]]);
					ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
				}
				else{
					if(S[s_ref[r]].flag==0){
						ln_b_fac -= l_dif * prm.U[sigma][s_ref[r]];
					}
				}
			}
		}
	}
	else{  // segment shrinks
		ln_b_fac = (-l_dif) * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
				if(S[s_ref[r]].k){
					double l_over = tau_overlap(tau1, F->tau1[i_tau1], S[s_ref[r]]);
					ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
				}
				else{
					if(S[s_ref[r]].flag==0){
						ln_b_fac += (-l_dif) * prm.U[sigma][s_ref[r]];
					}
				}
			}
		}
	}

	double b_fac = exp(ln_b_fac);

	//
	// determinant
	//
	double Delta_tau1[F->k];
	for(int i=0; i<F->k; i++){
		Delta_tau1[i] = G0_calc_interp(Delta[sigma], F->tau2[i], tau1);
	}

	// double lambda = shift1_lambda(*F, i_tau1, Delta_tau1);
	double lambda = F->D.shift1_lambda(i_tau1, Delta_tau1);


	//
	// phonon part
	//
	double fac_ph = 1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 0;

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau1[i] - tau1))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau1)))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau1[i_tau1]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - tau1))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau1)))
				           + exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau1[i_tau1]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau1[i_tau1])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= exp(-prm.w_ph * fabs(F->tau1[i_tau1] - tau1))
		           + exp(-prm.w_ph * (prm.beta - fabs(F->tau1[i_tau1] - tau1)));
		ln_fac_ph += 1.0 + exp(-prm.w_ph * prm.beta);

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = 0;

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau1)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau1[i_tau1])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau1)))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau1[i_tau1])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau1[i_tau1] - tau1)));
		ln_fac_ph += cosh(prm.w_ph * prm.beta*0.5);

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = lambda * b_fac;
	double prob = lambda * b_fac * fac_ph;
	if(flag_change==1)  prob = -prob;

	#if TEST_MODE
	printf("\nshift_tau1  (sigma=%d)\n", sigma);
	printf("  i_tau1=%d  tau1=%8.5lf to %8.5lf\n", i_tau1, F->tau1[i_tau1], tau1);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf("  prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_shift, ACCPT_sub.n_reject_shift) ){
		//
		// update tau-array
		//
		F->tau1[i_tau1] = tau1;

		//
		// update mat_M
		//
		// shift1_mat_M(*F, i_tau1, Delta_tau1, lambda);
		F->D.shift1_mat_M(i_tau1, Delta_tau1, lambda);

		//
		// rotate rows
		//
		if( flag_change==1 ){
			if( F->flag ){  // flag==1: tau1[i] > tau2[i],  rotate upward
				// rotate_upward(F->tau1, F->k);
				F->rotate_upward_tau1();
				// rotate_rows_upward(F->D.mat_M, F->k);
				F->D.rotate_rows_upward();
			}
			else{  // flag==1: tau1[i] > tau2[i],  rotate downward
				// rotate_downward(F->tau1, F->k);
				F->rotate_downward_tau1();
				// rotate_rows_downward(F->D.mat_M, F->k);
				F->D.rotate_rows_downward();
			}
		}

		F->tau1[F->k] = F->tau1[0] + prm.beta;

		F->flag = (F->flag + flag_change)%2;
		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}

static void shift_tau2_rand(int sigma)
{
	if(S[sigma].k)  shift_tau2(sigma, rand_int(S[sigma].k));;
}

static void shift_tau2_all(Operators *S, int sigma)
{
	if(S[sigma].flag){  // flag==1: tau1[i] > tau2[i],  execute in descending order
		for(int i=S[sigma].k-1; i>=0; i--){
			shift_tau2(sigma, i);
		}
	}
	else{  // flag==0: tau1[i] < tau2[i],  execute in ascending order
		for(int i=0; i<S[sigma].k; i++){
			shift_tau2(sigma, i);
		}
	}

	#if TEST_MODE
	printf("\n");
	for(int i=0; i<S[sigma].k; i++){
		printf("%3d  %8.5lf  %8.5lf\n", i, S[sigma].tau1[i], S[sigma].tau2[i]);
	}
	printf("  flag=%d\n", S[sigma].flag);
	#endif
}

//
// shift tau2 (f-creation op., c-annihilation op.)
//
static void shift_tau2(int sigma, int i_tau2)
{
	Operators *F = &S[sigma];

	int s_ref[N_S-1];
	for(int r=0; r<N_S-1; r++){
		if(r<sigma)  s_ref[r] = r;
		else         s_ref[r] = r+1;
	}

// 	if(F->k==0)  return;
//
// 	int i_tau2 = rand_int(F->k);

// 	double l_over, l_total;
// 	double b_fac;  // Boltzmann factor
	double l_max, l_dif, tau2;
	int flag_change=0;
// 	int flag_direction;  // 0: segment shrinks,  1: segment extends


	if( F->flag ){  // flag==1: tau1[i] > tau2[i]
		if(i_tau2)  l_max = F->tau1[i_tau2] - F->tau1[i_tau2-1];
		else  l_max = F->tau1[F->k] - F->tau1[F->k-1];
		if( prm.UINF ){
			for(int r=0; r<N_S-1; r++){
				// cond_op *p_S = &S[s_ref[r]];
				Operators *p_S = &S[s_ref[r]];
				if( p_S->k ){
					// int i_tau1_temp = tau_order(p_S->tau1, p_S->k, F->tau2[i_tau2]);
					int i_tau1_temp = tau_order(p_S->tau1, F->tau2[i_tau2]);
					i_tau1_temp = (i_tau1_temp - 1 + p_S->k) % p_S->k;
					double l_max_temp = F->tau1[i_tau2] - p_S->tau1[i_tau1_temp];
					if( l_max_temp < 0 )  l_max_temp += prm.beta;
					l_max = MIN(l_max, l_max_temp);
				}
			}
		}

		tau2 = F->tau1[i_tau2] - rand_tau_l(l_max);
		l_dif = tau2 - F->tau2[i_tau2];

		if( tau2 < 0 ){
			flag_change = 1;
			tau2 += prm.beta;
		}
	}
	else{  // flag==0: tau1[i] < tau2[i]
		l_max = F->tau1[i_tau2+1] - F->tau1[i_tau2];
		double tau1_temp = F->tau1[i_tau2];
		if( prm.UINF ){
			for(int r=0; r<N_S-1; r++){
				// cond_op *p_S = &S[s_ref[r]];
				Operators *p_S = &S[s_ref[r]];
				if( p_S->k ){
					// int i_tau1_temp = tau_order(p_S->tau1, p_S->k, F->tau2[i_tau2]);
					int i_tau1_temp = tau_order(p_S->tau1, F->tau2[i_tau2]);
					i_tau1_temp = (i_tau1_temp - 1 + p_S->k) % p_S->k;
					double l_max_temp = F->tau1[i_tau2+1] - p_S->tau1[i_tau1_temp];
					if( l_max_temp < 0 )  l_max_temp += prm.beta;
					if( l_max > l_max_temp ){
						l_max = l_max_temp;
						tau1_temp = p_S->tau1[i_tau1_temp];
					}
				}
			}
		}

// 		tau2 = F->tau1[i_tau2] + rand_tau_l(l_max);
		tau2 = tau1_temp + rand_tau_l(l_max);
		l_dif = tau2 - F->tau2[i_tau2];

		if( tau2 > prm.beta ){
			flag_change = 1;
			tau2 -= prm.beta;
		}
	}

	double ln_b_fac;

	if( l_dif < 0 ){ // segment extends
		ln_b_fac = -(-l_dif) * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
				if(S[s_ref[r]].k){
					double l_over = tau_overlap(tau2, F->tau2[i_tau2], S[s_ref[r]]);
					ln_b_fac -= l_over * prm.U[sigma][s_ref[r]];
				}
				else{
					if(S[s_ref[r]].flag==0){
						ln_b_fac -= (-l_dif) * prm.U[sigma][s_ref[r]];
					}
				}
			}
		}
	}
	else{  // segment shrinks
		ln_b_fac = l_dif * prm.ef[sigma];
		if( prm.UINF==0 ){
			for(int r=0; r<N_S-1; r++){
				if(S[s_ref[r]].k){
					double l_over = tau_overlap(F->tau2[i_tau2], tau2, S[s_ref[r]]);
					ln_b_fac += l_over * prm.U[sigma][s_ref[r]];
				}
				else{
					if(S[s_ref[r]].flag==0){
						ln_b_fac += l_dif * prm.U[sigma][s_ref[r]];
					}
				}
			}
		}
	}

	double b_fac = exp(ln_b_fac);

	//
	// determinant
	//
	double Delta_tau2[F->k];
	for(int i=0; i<F->k; i++){
// 		Delta_tau2[i] = calc_Delta_interp(tau2, F->tau1[i]);
		Delta_tau2[i] = G0_calc_interp(Delta[sigma], tau2, F->tau1[i]);
	}

	// double lambda = shift2_lambda(*F, i_tau2, Delta_tau2);
	double lambda = F->D.shift2_lambda(i_tau2, Delta_tau2);


	//
	// phonon part
	//
	double fac_ph = 1.0;

	#if PHONON
	if(prm.beta * prm.w_ph > 1){  // low-T

		double ln_fac_ph = 0;

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += exp(-prm.w_ph * fabs(S[s].tau2[i] - tau2))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - tau2)))
				           - exp(-prm.w_ph * fabs(S[s].tau2[i] - F->tau2[i_tau2]))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau2[i] - F->tau2[i_tau2])))
				           - exp(-prm.w_ph * fabs(S[s].tau1[i] - tau2))
				           - exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - tau2)))
				           + exp(-prm.w_ph * fabs(S[s].tau1[i] - F->tau2[i_tau2]))
				           + exp(-prm.w_ph * (prm.beta - fabs(S[s].tau1[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= exp(-prm.w_ph * fabs(F->tau2[i_tau2] - tau2))
		           + exp(-prm.w_ph * (prm.beta - fabs(F->tau2[i_tau2] - tau2)));
		ln_fac_ph += 1.0 + exp(-prm.w_ph * prm.beta);

		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * (1.0 - exp(-prm.beta * prm.w_ph)));

		fac_ph = exp(ln_fac_ph);
	}
	else{  // high-T

		double ln_fac_ph = 0;

		for(int s=0; s<N_S; s++){
			for(int i=0; i<S[s].k; i++){
				ln_fac_ph += cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - tau2)))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau2[i] - F->tau2[i_tau2])))
				           - cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - tau2)))
				           + cosh(prm.w_ph * (prm.beta*0.5 - fabs(S[s].tau1[i] - F->tau2[i_tau2])));
			}
		}

		// subtract unnecessary contribution
		ln_fac_ph -= cosh(prm.w_ph * (prm.beta*0.5 - fabs(F->tau2[i_tau2] - tau2)));
		ln_fac_ph += cosh(prm.w_ph * prm.beta*0.5);

	// 	double fac_ph = prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));
		ln_fac_ph *= -prm.g * prm.g / (prm.w_ph * prm.w_ph * sinh(prm.beta * prm.w_ph * 0.5));

		fac_ph = exp(ln_fac_ph);
	}
	#endif // PHONON


// 	double prob = lambda * b_fac;
	double prob = lambda * b_fac * fac_ph;
	if(flag_change==1)  prob = -prob;

	#if TEST_MODE
	printf("\nshift_tau2  (sigma=%d)\n", sigma);
	printf("  i_tau2=%d  tau2=%8.5lf to %8.5lf\n", i_tau2, F->tau2[i_tau2], tau2);
	printf("  l_max=%8.5lf\n", l_max);
	printf("  lambda=%.4lf\n", lambda);
	printf("  b_fac=%.3e\n", b_fac);
	printf(" prob=%.5lf\n", prob);
	#endif


	if( metropolis_abs(fabs(prob), ACCPT_sub.n_accept_shift, ACCPT_sub.n_reject_shift) ){
		//
		// update tau-array
		//
		F->tau2[i_tau2] = tau2;

		//
		// update mat_M
		//
		// shift2_mat_M(*F, i_tau2, Delta_tau2, lambda);
		F->D.shift2_mat_M(i_tau2, Delta_tau2, lambda);

		//
		// rotate columns
		//
		if( flag_change==1 ){
			if( F->flag ){  // flag==1: tau1[i] > tau2[i],  rotate downward
				// rotate_downward(F->tau2, F->k);
				F->rotate_downward_tau2();
				// rotate_columns_downward(F->D.mat_M, F->k);
				F->D.rotate_columns_downward();
			}
			else{  // flag==1: tau1[i] > tau2[i],  rotate upward
				// rotate_upward(F->tau2, F->k);
				F->rotate_upward_tau2();
				// rotate_columns_upward(F->D.mat_M, F->k);
				F->D.rotate_columns_upward();
			}
		}

		F->tau2[F->k] = F->tau2[0] + prm.beta;

		F->flag = (F->flag + flag_change)%2;
		if(prob<0)  w_sign = -w_sign;

		#if TEST_MODE
		printf("\n");
		for(int i=0; i<F->k; i++){
			printf("%3d  %8.5lf  %8.5lf\n", i, F->tau1[i], F->tau2[i]);
		}
		printf("  flag=%d\n", F->flag);
		#endif
	}
}


// This function works for all values of k including k=0
static double tau_length(Operators &F)
{
	double l=0;

	for(int i=0; i<F.k; i++){
		l += F.tau1[i] - F.tau2[i];
	}
	if(F.flag==0)  l += prm.beta;

	return(l);
}

// static int tau_overlap(double tau_from, double tau_to, Operators &F, double *l_over, double *l_total)
static double tau_overlap(double tau_from, double tau_to, Operators &F)
{
	if( F.k == 0 ){
		double l_over = 0;
		if( F.flag == 0 ){  // |1>
			l_over = tau_to - tau_from;
			if( l_over < 0 )  l_over += prm.beta;
		}
		return l_over;
	}

	// double *tau1, *tau2;
	std::vector<double> *tau1, *tau2;
	double *l1, *l2, l_comp, l_over, l_total;
	int n=0;

	if(F.flag){
		tau1 = &F.tau1;
		tau2 = &F.tau2;
		l1 = &l_over;
		l2 = &l_comp;
	}
	else{
		tau1 = &F.tau2;
		tau2 = &F.tau1;
		l1 = &l_comp;
		l2 = &l_over;
	}

	*l1 = 0;

	// int i_from=tau_order(tau2, F.k, tau_from);
	int i_from=tau_order(*tau2, tau_from);
	// int i_to=tau_order(tau2, F.k, tau_to);
	int i_to=tau_order(*tau2, tau_to);

	#if TEST_MODE==2
	printf("\ntau_overlap\n");
	printf("  i_from=%d, i_to=%d\n", i_from, i_to);
	#endif

	if(tau_from < tau_to){

		for(int i=i_from; i<i_to; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if TEST_MODE==2
			printf(" (*tau1)[%d]-tau2[%d]=%.5lf\n", i, i, (*tau1)[i] - tau2[i]);
			#endif
		}

		if(i_from){  // !=0
			double under_estim = (*tau1)[i_from-1] - tau_from;
			if( under_estim > 0 ){
				*l1 += under_estim;
				n += 1;

				#if TEST_MODE==2
				printf(" under_estim=%.5lf\n", under_estim);
				#endif
			}
		}

		if(i_to){  // !=0
			double over_estim = (*tau1)[i_to-1] - tau_to;
			if( over_estim > 0 ){
				*l1 -= over_estim;
				n -= 1;

				#if TEST_MODE==2
				printf(" over_estim=%.5lf\n", over_estim);
				#endif
			}
		}

// 		if(F.flag==0)  *l1 = tau_to - tau_from - *l1;
		l_total = tau_to - tau_from;
		*l2 = l_total - *l1;
	}
	else{

		for(int i=i_from; i<F.k; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if TEST_MODE==2
			printf(" tau1[%d]-tau2[%d]=%.5lf\n", i, i, (*tau1)[i] - (*tau2)[i]);
			#endif
		}

		if(i_from){  // !=0
			double under_estim = (*tau1)[i_from-1] - tau_from;
			if( under_estim > 0 ){
				*l1 += under_estim;
				n += 1;

				#if TEST_MODE==2
				printf(" under_estim=%.5lf\n", under_estim);
				#endif
			}
		}

		for(int i=0; i<i_to; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if TEST_MODE==2
			printf(" (*tau1)[%d]-(*tau2)[%d]=%.5lf\n", i, i, (*tau1)[i] - (*tau2)[i]);
			#endif
		}

		if(i_to){  // !=0
			double over_estim = (*tau1)[i_to-1] - tau_to;
			if( over_estim > 0 ){
				*l1 -= over_estim;
				n -= 1;

				#if TEST_MODE==2
				printf(" over_estim=%.5lf\n", over_estim);
				#endif
			}
		}

// 		if(F.flag==0)  *l1 = prm.beta + tau_to - tau_from - *l1;
		l_total = prm.beta + tau_to - tau_from;
		*l2 = l_total - *l1;
	}

	#if TEST_MODE==2
	printf("  l_overlap=%.5lf  l_total=%.5lf\n", l_over, l_total);
// 	if(F.flag)  printf("  l_overlap=%.5lf  l_total=%.5lf\n", *l1, *l_total);
// 	else  printf("  l_overlap=%.5lf  l_total=%.5lf\n", *l2, *l_total);
	printf("  n_overlap=%d\n", n);
	#endif

// 	return(n);
	return( l_over );
}

//============================================================================
