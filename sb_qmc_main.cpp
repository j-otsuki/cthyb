static const char tag[64] = "v1.01";

#include "ct_qmc_share.h"
#include "sb_qmc.h"
#include "green_func_0.h"
#include "pade.h"
#include <gsl/gsl_integration.h>
// #include "g_noninteract.h"


struct sb_qmc_params prm;
// double prm_ave_n, prm_chem_pot;
// double prm_V_sqr[N_S], prm_V_sqr_imp[N_S], prm_D, prm_f_disp;
struct num_mc n_mc;

// struct statistics *STAT;
struct phys_quant *PQ;
struct single_particle *SP;
struct two_particle *TP;
struct vertex *VX;

static complex<double> Delta_omega[N_TAU/2];
static complex<double> D0_omega[N_TAU/2+1];
static double prm_V_sqr;
static double prm_D;
static double prm_gsqr_ch, prm_gsqr_z, prm_gsqr_xy;
static double prm_ef;
static double prm_h;

static double PQ_energy_tot, PQ_energy_tot_err, PQ_energy_kin, PQ_energy_int;

static int prm_dos_boson;
// 0: DOS = delta(w-w0);
// 1: DOS \propto w^{gamma}  (with cutoff `w_cutoff')
static double prm_w0;  // if dos_boson == 0
static double prm_w_cutoff, prm_gamma;  // if dos_boson == 1

static int N_PADE_SP = MIN(N_TAU/2, 8192);
// Number of frequencies used for Pade approximation of Gf (N_PADE_SP^2 memory)
// MAX: N_TAU/2 (all data),  MIN: 0 (do not evaluate)

static int N_X, X_MESH;
static double X_MIN, X_MAX;

static int my_rank=0, process_num=1;

int i_program;
char str_program[][128]={
	"Parameters fixed",
	"T varies",
	"V^2 varies",
	"g^2 (= g_xy = g_z) varies",
	"g_xy^2 varies",
	"g_z^2 varies",
	"ef varies",
	"h varies",
	"w0 varies",
	"gamma varies",
};

int flag_tp;
char str_flag[][128]={
	"Green function",  // 0
	"Green function, susceptibility",  // 1
	"Green function, susceptibility, vertex",  // 2
	"free-energy by Wang-Landau sampling", // 3
};


int main_sub(int num, double x);
void main_x();
void init_params();

void print_pq(int num, double x);
void print_thermo(int, double x);
void print_energy(int);
void print_occup();
void print_Gf0(int num);
void print_statistics(int);
void print_single_particle(int);
void print_two_particle_pm(int);
void print_two_particle_zz(int);
void print_D0(int);
void print_vertex(int);

static void eval_energy(double beta, double cutoff, double gamma);



main(int argc, char* argv[])
{
	#if QMC_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	#endif // QMC_MPI
	
	if(my_rank==0)  printf("SB_QMC_MAIN  %s\n", tag);
	
	init_params();
	
	sbqmc_init(&PQ, &SP, &TP, &VX);
// 	sbqmc_init();
	sbqmc_set_nmc(n_mc);
	
// 	G = (dmft_green_func *)malloc(N_S*sizeof(dmft_green_func));
// 	for(int s=0; s<N_S; s++){
// 		G[s].Gf_perio = G[s].tmat_perio;
// 		G[s].F = G[s].self_c;
// 	}
	
	
	if( i_program == 0 )  main_sub(0, 0);
	else                  main_x();
	
// 	free(G);
	
// 	delete [] G0_omega;
// 	delete [] store_Gc_cav_omega;
	
	sbqmc_final();
	
	#if QMC_MPI
	MPI_Finalize();
	#endif // QMC_MPI
}

// return 0: computed
// return 1: not computed ( ave_sign < MIN_AVE_SIGN )
int main_sub(int num, double x)
{
	clock_t time_start = clock();
	
	// init G0
	int flag_V = 1;
	if( prm_V_sqr == 0 )  flag_V = 0;
	
	if( flag_V ){
		G0_omega_calc(Delta_omega, prm.beta, prm_D);
		for(int i=0; i<N_TAU/2; i++){
			Delta_omega[i] *= prm_V_sqr;
		}
	// 	G0_omega_calc(G0_omega, N_S, prm.beta, prm_D, flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
		if(my_rank==0){
// 			print_Gf0(num);
		}
	}
	
	// init D0
	int flag_g = 1;
	if( prm_gsqr_xy == 0 )  flag_g = 0;
	
	if( prm_dos_boson == 0 ){
		D0_omega_calc0(D0_omega, prm.beta, prm_w0);
	}
	else if( prm_dos_boson == 1 ){
		D0_omega_calc1(D0_omega, prm.beta, prm_w_cutoff, prm_gamma);
	}
	if(my_rank==0){
// 		print_D0(num);
	}
	
	complex<double> D0_ch_omega[N_TAU/2+1], D0_sp_omega[N_TAU/2+1], D0_pm_omega[N_TAU/2+1];
	for(int i=0; i<=N_TAU/2; i++){
		D0_ch_omega[i] = D0_omega[i] * prm_gsqr_ch;
		D0_sp_omega[i] = D0_omega[i] * prm_gsqr_z;
		D0_pm_omega[i] = D0_omega[i] * prm_gsqr_xy;
	}
	complex<double> *p_D0_omega[3] = {D0_ch_omega, D0_sp_omega, D0_pm_omega};
	if( prm_gsqr_ch == 0 )  p_D0_omega[0] = NULL;
	if( prm_gsqr_z  == 0 )  p_D0_omega[1] = NULL;
	if( prm_gsqr_xy == 0 )  p_D0_omega[2] = NULL;
	
	complex<double> *p_Delta[2] = {Delta_omega, Delta_omega};
	double v_sqr[2] = {prm_V_sqr, prm_V_sqr};
// 	double fac_int[3] = {1, 1, 1};  // Pauli
	double fac_int[3] = {1, 0.25, 0.25};  // spin
	double fac_chi[3] = {1, 1, 1};  // Pauli
// 	double fac_chi[3] = {1, 0.25, 0.25};  // spin
// 	sbqmc_set_input1(prm, p_Delta, v_sqr, p_D0_omega);  // spin
// 	sbqmc_set_input2(prm, p_Delta, v_sqr, p_D0_omega);  // sigma
	sbqmc_set_input(prm, p_Delta, v_sqr, p_D0_omega, fac_int, fac_chi);
	
	if( flag_tp < 3 ){
		int info = sbqmc_eval(flag_tp);
		if( info ){
			printf("\n sbqmc_eval() returned %d\n", info);
			exit(info);
		}
		
		if(my_rank==0){
			print_pq(num, x);
			print_statistics(num);
			print_occup();
			
			if( flag_V )  print_single_particle(num);
			if( flag_g )  print_two_particle_pm(num);
			
			if(flag_tp >= 1){
				print_two_particle_zz(num);
			}
			if(flag_tp >= 2){
				print_vertex(num);
			}
			
// 			if( LIMIT_PARAM == 2 ){
			if( ! flag_V ){  // spin system
				if( prm_dos_boson == 1 ){
					eval_energy(prm.beta, prm_w_cutoff, prm_gamma);
					print_energy(num);
				}
			}
		}
	}
	else{
		int info = sbqmc_eval_thermo();
		
		if(my_rank==0){
			print_thermo(num, x);
			print_statistics(num);
		}
	}
	
	if(my_rank==0){
		clock_t time_end = clock();
// 		print_time(time_start, time_end);
// 	// 	print_time(time_start, time_end, LOG_FILE);
	}
}

void main_x()
{
	clock_t time_start = clock();
	
	double x[N_X];
// 	double T_energy[N_X], T_energy_err[N_X];
	
	mesh_init(X_MESH, x, X_MIN, X_MAX, N_X, my_rank);
	
	
	for(int ix=0; ix<N_X; ix++){
		
		if(my_rank==0){
			char str_t[128];
			sprintf(str_t, "\n===   ix=%d  x=%lf\n", ix, x[ix]);
			printf("%s", str_t);
			sbqmc_fprint_log(str_t);
		}
		
		switch(i_program){
			case 1:
				prm.beta = 1.0 / x[ix];
				break;
			case 2:
				prm_V_sqr = x[ix];
				break;
			case 3:
				prm_gsqr_xy = prm_gsqr_z = x[ix];
				break;
			case 4:
				prm_gsqr_xy = x[ix];
				break;
			case 5:
				prm_gsqr_z = x[ix];
				break;
			case 6:
				prm_ef = x[ix];
				break;
			case 7:
				prm_h = x[ix];
				break;
			case 8:
				prm_w0 = x[ix];
				break;
			case 9:
				prm_gamma = x[ix];
				break;
			default:
				exit(0);
		}
		prm.ef[0] = prm_ef - prm_h / 2.;
		prm.ef[1] = prm_ef + prm_h / 2.;
		
		main_sub(ix, x[ix]);
		
// 		struct phys_quant{
// 			double ave_sign, ave_sign_err;
// 			double Z_V[2][N_K], Z_V_err[2][N_K];  // V
// 			double Z_g[N_K], Z_g_err[N_K];  // g
// 			double ave_Z_V[2], ave_Z_V_err[2];
// 			double ave_Z_g, ave_Z_g_err;
// 			double occup_tot, occup_tot_err;
// 			double occup_dif, occup_dif_err;
// 			double occup_ud[2], occup_ud_err[2];
// 			double stat_suscep_sp, stat_suscep_sp_err;
// 			double stat_suscep_ch, stat_suscep_ch_err;
// 			double stat_suscep_pm, stat_suscep_pm_err;
// 		};
		if(my_rank==0){
		}
	}
	
	
	if(my_rank==0){
		printf("\n===\n");
		
		clock_t time_end = clock();
		print_time(time_start, time_end);
	// 	print_time(time_start, time_end, LOG_FILE);
	}
	
}


void init_params()
{
	FILE *fp;
	if( (fp = fopen("sb_qmc_param.init", "r")) == NULL ){
		if(my_rank==0){
			printf("'sb_qmc_param.init' not found\n");
		}
		exit(0);
	}
	
	double data[16];
	
	if( read_data(fp, data) != 1 )  exit(0);
	i_program = (int)data[0];
	
	if( read_data(fp, data) != 1 )  exit(0);
	flag_tp = (int)data[0];
	
	
	if( read_data(fp, data) != 1 )  exit(0);
	prm_D = (double)data[0];
	
	if( read_data(fp, data) != 1 )  exit(0);
	prm.beta = 1./ (double)data[0];
	
	if( read_data(fp, data) != 2 )  exit(0);
	prm_ef = (double)data[0];
	prm_h = (double)data[1];
	prm.ef[0] = prm_ef - prm_h / 2.;
	prm.ef[1] = prm_ef + prm_h / 2.;
	
	if( read_data(fp, data) != 1 )  exit(0);
	prm.U = (double)data[0];
	
	
	if( read_data(fp, data) != 1 )  exit(0);
	prm_V_sqr = (double)data[0];
	
	if( read_data(fp, data) != 1 )  exit(0);
	prm_gsqr_ch = (double)data[0];
	
	switch( read_data(fp, data) ){
		case 1:
			prm_gsqr_z = prm_gsqr_xy = (double)data[0];
			break;
		case 2:
			prm_gsqr_xy = (double)data[0];
			prm_gsqr_z = (double)data[1];
			break;
		default:
			exit(0);
	}
	
	int n_data;
	n_data = read_data(fp, data);
	prm_dos_boson = (int)data[0];
	if( prm_dos_boson == 0 ){
		if(n_data != 2)  exit(0);
		prm_w0 = (double)data[1];
	}
	else if( prm_dos_boson == 1 ){
		if(n_data != 3)  exit(0);
		prm_w_cutoff = (double)data[1];
		prm_gamma = (double)data[2];
	}
	else{
		exit(0);
	}
	
	
	if(my_rank==0){
		printf("\n %d: %s\n", i_program, str_program[i_program]);
		printf(" %d: %s\n", flag_tp, str_flag[flag_tp]);
		
		printf("\n DOS = %d\n", DOS);
		printf(" D = %.4e\n", prm_D);
		printf(" T = %.4e  (beta = %.4e)\n", 1.0 / prm.beta, prm.beta);
		printf(" V^2 = %.4e\n", prm_V_sqr);
		printf(" g_ch = %.4e\n", prm_gsqr_ch);
		printf(" g_xy = %.4e  g_z = %.4e\n", prm_gsqr_xy, prm_gsqr_z);
		printf(" ef = %.4e  h = %.4e\n", prm_ef, prm_h);
		if( UINF )  printf(" U = infinity\n");
		else  printf(" U = %.4e\n", prm.U);
		printf(" dos_boson = %d", prm_dos_boson);
		if( prm_dos_boson == 0 ){
			printf("  w0 = %.4e\n", prm_w0);
		} else {
			printf("  w_cutoff = %.4e  gamma = %.4e\n", prm_w_cutoff, prm_gamma);
		}
	}
	
	if( read_data(fp, data) != 2 )  exit(0);
	n_mc.bin = (int)data[0];
	n_mc.msr = (int)data[1];
	
	if( read_data(fp, data) != 5 )  exit(0);
	n_mc.seg = (int)data[0];
	n_mc.spin = (int)data[1];
	n_mc.boson = (int)data[2];
	n_mc.shift = (int)data[3];
	n_mc.segsp = (int)data[4];
	
	if(my_rank==0){
		printf("\n Numbers of sampling/updates\n");
		printf("  bin = %d  msr = %d\n", n_mc.bin, n_mc.msr);
		printf("  seg = %d  spin = %d", n_mc.seg, n_mc.spin);
		printf("  boson = %d  shift = %d", n_mc.boson, n_mc.shift);
		printf("  segsp = %d\n", n_mc.segsp);
	}
	
	if(i_program){
		if( read_data(fp, data) != 4 )  exit(0);
		X_MESH = (int)data[0];
		N_X = (int)data[1];
		X_MIN = data[2];
		X_MAX = data[3];
	}
	
	fclose(fp);
	
	fflush(stdout);
}


void print_pq(int num, double x)
{
	FILE *fp;
	
	char filename[128];
	sprintf(filename, DATA_DIR "xx.dat");
	
	if(num==0)  fp=fopen(filename, "w");
	else        fp=fopen(filename, "a");
	
	// 1-4
	fprintf(fp, "%.5e", x);
	fprintf(fp, " %.6e", 1./ prm.beta);
	fprintf(fp, " %.6e %.6e", PQ->ave_sign, PQ->ave_sign_err);
	// 5-18
	for(int s=0; s<2; s++){
		fprintf(fp, " %.6e %.6e", PQ->ave_Z_V[s], PQ->ave_Z_V_err[s]);
	}
	fprintf(fp, " %.6e %.6e", PQ->ave_Z_g, PQ->ave_Z_g_err);
	fprintf(fp, " %.6e %.6e", PQ->occup_tot, PQ->occup_tot_err);
	fprintf(fp, " %.6e %.6e", PQ->occup_dif, PQ->occup_dif_err);
// 	for(int s=0; s<2; s++){
// 		fprintf(fp, " %.6e %.6e", PQ->occup_ud[s], PQ->occup_ud_err[s]);
// 	}
	fprintf(fp, " %.6e %.6e", PQ->occup_dbl, PQ->occup_dbl_err);
	fprintf(fp, " ? ?");
	// 19-24
	fprintf(fp, " %.6e %.6e", PQ->stat_suscep_sp, PQ->stat_suscep_sp_err);
	fprintf(fp, " %.6e %.6e", PQ->stat_suscep_ch, PQ->stat_suscep_ch_err);
	fprintf(fp, " %.6e %.6e", PQ->stat_suscep_pm, PQ->stat_suscep_pm_err);
	// 25-48
	for(int s=0; s<2; s++){
		fprintf(fp, " %.6e %.6e", real(SP[s].Gf_omega[0]), imag(SP[s].Gf_omega[0]));
		fprintf(fp, " %.6e %.6e", real(SP[s].Gf_omega[1]), imag(SP[s].Gf_omega[1]));
	}
	for(int s=0; s<2; s++){
		fprintf(fp, " %.6e %.6e", real(SP[s].self_f_1[0]), imag(SP[s].self_f_1[0]));
		fprintf(fp, " %.6e %.6e", real(SP[s].self_f_1[1]), imag(SP[s].self_f_1[1]));
	}
	if( UINF ){
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_f_2[0]), imag(SP[s].self_f_2[0]));
			fprintf(fp, " %.6e %.6e", real(SP[s].self_f_2[1]), imag(SP[s].self_f_2[1]));
		}
	}
	else{
		fprintf(fp, " ? ? ? ? ? ? ? ?");
	}
	// 49-52
	// chi(beta/2)
	fprintf(fp, " %.6e %.6e", TP->chi_sp_tau[N_TP], TP->chi_sp_tau_err[N_TP]);
	fprintf(fp, " %.6e %.6e", TP->chi_pm_tau[N_TAU/2], TP->chi_pm_tau_err[N_TAU/2]);
	//53-56
	fprintf(fp, " %.6e %.6e", PQ->m_z_sqr, PQ->m_z_sqr_err);
	fprintf(fp, " %.6e %.6e", PQ->m_z_pw4, PQ->m_z_pw4_err);
	//57-60
	fprintf(fp, " %.6e %.6e", PQ->m_xy_sqr, PQ->m_xy_sqr_err);
	fprintf(fp, " %.6e %.6e", PQ->m_xy_pw4, PQ->m_xy_pw4_err);
	//61-64
	fprintf(fp, " %.6e %.6e", PQ->m_rt_sqr, PQ->m_rt_sqr_err);
	fprintf(fp, " %.6e %.6e", PQ->m_rt_pw4, PQ->m_rt_pw4_err);
	
	// 65-68
	for(int i=0; i<4; i++)  fprintf(fp, " %.6e", real(TP->chi_sp_omega[i]));
	// 69-72
	for(int i=0; i<4; i++)  fprintf(fp, " %.6e", real(TP->chi_pm_omega[i]));
	// 73-74
	fprintf(fp, " %.6e %.6e", PQ->eff_mom, PQ->stat_suscep_reg);
	
	fprintf(fp, "\n");
	fclose(fp);
	printf("\n '%s' updated\n", filename);
}

void print_thermo(int num, double x)
{
	printf("\n free_energy = %.5e +- %.4e\n", PQ->free_energy, PQ->free_energy_err);
	
	FILE *fp;
	char filename[128];
	sprintf(filename, DATA_DIR "xx-thermo.dat");
	
	if(num==0)  fp=fopen(filename, "w");
	else        fp=fopen(filename, "a");
	
// 	fprintf(fp, "%.5e", x);
	fprintf(fp, "%.8e", 1./ prm.beta);
	fprintf(fp, " %.8e %.8e", PQ->free_energy, PQ->free_energy_err);
	fprintf(fp, "\n");
	fclose(fp);
	printf("\n '%s' updated\n", filename);
}

void print_energy(int num)
{
	FILE *fp;
	
	char filename[128];
	sprintf(filename, DATA_DIR "xx-energy.dat");
	
	if(num==0)  fp=fopen(filename, "w");
	else        fp=fopen(filename, "a");
	
// 	fprintf(fp, "%.5e", x);
	fprintf(fp, "%.8e", 1./ prm.beta);
	fprintf(fp, " %.8e %.8e", PQ_energy_tot, PQ_energy_tot_err);
	fprintf(fp, " %.8e %.8e", PQ_energy_kin, PQ_energy_int);
	fprintf(fp, "\n");
	fclose(fp);
	printf("\n '%s' updated\n", filename);
}

void print_statistics(int num)
{
	FILE *fp;
	char filename[128];
	
	sprintf(filename, DATA_DIR "%02d-stat.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_K; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<2; s++){
			if( PQ->Z_V[s][i] )
				fprintf(fp, " %.6e %.5e", PQ->Z_V[s][i], PQ->Z_V_err[s][i]);
			else
				fprintf(fp, " ? ?");
		}
		if( PQ->Z_g[i] )
			fprintf(fp, " %.6e %.5e", PQ->Z_g[i], PQ->Z_g_err[i]);
		else
			fprintf(fp, " ? ?");
		fprintf(fp, "\n");
	}
	printf("\n '%s'\n", filename);
	fclose(fp);
}

void print_occup()
{
	printf("\n occupation\n");
	
	for(int s=0; s<2; s++){
		printf("  %3d: %lf +- %lf\n", s, PQ->occup_ud[s], PQ->occup_ud_err[s]);
	}
	printf("  tot: %lf +- %lf\n", PQ->occup_tot, PQ->occup_tot_err);
	printf("  dif: %lf +- %lf\n", PQ->occup_dif, PQ->occup_dif_err);
	printf("  dbl: %lf +- %lf\n", PQ->occup_dbl, PQ->occup_dbl_err);
	
	printf("\n moment\n");
	printf("  <m^2> : %lf +- %lf\n", PQ->m_z_sqr, PQ->m_z_sqr_err);
	printf("  <m^4> : %lf +- %lf\n", PQ->m_z_pw4, PQ->m_z_pw4_err);
}

void print_single_particle(int num)
{
	FILE *fp;
	char filename[128];
	
	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}
	
	
	sprintf(filename, DATA_DIR "%02d-Gf_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.5e", tau[i]);
		
// 		for(int s=0; s<1; s++){
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", SP[s].Gf_tau[i], SP[s].Gf_tau_err[i]);
		}
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e", SP[s].G_self[i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-Gf_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].Gf_omega[i]), imag(SP[s].Gf_omega[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	complex<double> i_omega_f[N_PADE_SP];
	
	if(N_PADE_SP){
		complex<double> Gf_pade[2][N_PADE_SP];
		
		for(int i=0; i<N_PADE_SP; i++){
			double omega_f = (double)(2*i+1) * M_PI / prm.beta;
			i_omega_f[i] = IMAG * omega_f;
		}
		
		for(int s=0; s<2; s++){
			pade_init(i_omega_f, SP[s].Gf_omega, Gf_pade[s], N_PADE_SP);
		}
		
		sprintf(filename, DATA_DIR "%02d-Gf_pade.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<1000; i++){
			
			complex<double> w = -2 * prm_D + 4 * prm_D / 999. * (double)i;
			
			fprintf(fp, "%.6e", real(w));
			for(int s=0; s<2; s++){
				complex<double> temp_Gf = pade_calc(i_omega_f, Gf_pade[s], w, N_PADE_SP);
				fprintf(fp, " %.6e %.6e", real(temp_Gf), imag(temp_Gf));
			}
			
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
	
	
	sprintf(filename, DATA_DIR "%02d-self.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_f_1[i]), imag(SP[s].self_f_1[i]));
		}
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_f[i]), imag(SP[s].self_f[i]));
		}
		if( UINF ){
			for(int s=0; s<2; s++){
				fprintf(fp, " %.6e %.6e", real(SP[s].self_f_2[i]), imag(SP[s].self_f_2[i]));
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
}

void print_two_particle_pm(int num)
{
	FILE *fp;
	char filename[40];
	
	printf("\n static susceptibility\n");
	printf("  spin(+-) : %.5e +- %.4e\n", PQ->stat_suscep_pm, PQ->stat_suscep_pm_err);
	printf("           : %.5e (test)\n", real(TP->chi_pm_omega[0]));
	
	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}
	
	sprintf(filename, DATA_DIR "%02d-chi_pm_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.5e", tau[i]);
		fprintf(fp, " %.6e %.6e", TP->chi_pm_tau[i], TP->chi_pm_tau_err[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-chi_pm_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU/2; i++){
		fprintf(fp, "%d", i);
		fprintf(fp, " %.5e", real(TP->chi_pm_omega[i]));
		fprintf(fp, " %.5e", imag(TP->chi_pm_omega[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}

void print_two_particle_zz(int num)
{
	printf("\n static susceptibility\n");
	printf("  spin(zz) : %.5e +- %.4e\n", PQ->stat_suscep_sp, PQ->stat_suscep_sp_err);
	printf("  charge   : %.5e +- %.4e\n", PQ->stat_suscep_ch, PQ->stat_suscep_ch_err);
// 	printf("           : %.5e\n", real(TP->chi_sp_omega[0]));
	
	FILE *fp;
	char filename[40];
	
	sprintf(filename, DATA_DIR "%02d-chi_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TP; i++){
// 		double tau = (double)i * prm.beta / (double)(2*N_TP);
		fprintf(fp, "%.5e", TP->tau[i]);
		fprintf(fp, " %.5e %.4e", TP->chi_sp_tau[i], TP->chi_sp_tau_err[i]);
		fprintf(fp, " %.5e %.4e", TP->chi_ch_tau[i], TP->chi_ch_tau_err[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-chi_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TP2; i++){
		fprintf(fp, "%d", i);
		fprintf(fp, " %.5e", real(TP->chi_sp_omega[i]));
		fprintf(fp, " %.5e", real(TP->chi_ch_omega[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
/*	
	complex<double> chi_sp_pade[N_TP2], chi_ch_pade[N_TP2], i_omega_b[N_TP2];
	
	for(int i=0; i<N_TP2; i++)  i_omega_b[i] = IMAG * (double)(2*i) * M_PI / prm.beta;
	
	pade_init(i_omega_b, TP_sp->chi_omega, chi_sp_pade, N_TP2);
	pade_init(i_omega_b, TP_ch->chi_omega, chi_ch_pade, N_TP2);
	
	sprintf(filename, DATA_DIR "%02d-chi_pade.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<1000; i++){
		complex<double> w = 1.0 * prm_D / 999. * (double)i;
		
		complex<double> temp_chi_sp = pade_calc(i_omega_b, chi_sp_pade, w, N_TP2);
		complex<double> temp_chi_ch = pade_calc(i_omega_b, chi_ch_pade, w, N_TP2);
		
		fprintf(fp, "%.6e %.6e %.6e %.6e %.6e\n", real(w),
		 real(temp_chi_sp), imag(temp_chi_sp), real(temp_chi_ch), imag(temp_chi_ch));
	}
	fclose(fp);
	printf(" '%s'\n", filename);
*/	
	
	//
	// transverse susceptibility
	//
	#if CHI_TR
	
	sprintf(filename, DATA_DIR "%02d-chi_tr_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=2*N_TP2; i++){
		fprintf(fp, "%.5e", (double)i * prm.beta / (double)(2*N_TP2));
		fprintf(fp, " %.5e %.4e", TP_tr->chi_tau1[i], TP_tr->chi_tau1_err[i]);
		fprintf(fp, " %.5e %.4e", TP_tr->chi_tau2[i], TP_tr->chi_tau2_err[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-chi_tr_omega-t.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TP2; i++){
		fprintf(fp, "%d", i);
		fprintf(fp, " %.5e %.5e", real(TP_tr->chi_omega1[i]), real(TP_tr->chi_omega1[i]));
		fprintf(fp, " %.5e %.5e", real(TP_tr->chi_omega2[i]), real(TP_tr->chi_omega2[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	sprintf(filename, DATA_DIR "%02d-chi_tr_omega-w.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TP_W; i++){
		fprintf(fp, "%d", i);
		fprintf(fp, " %.5e %.5e", real(TP_tr->chi_omega[i]), real(TP_tr->chi_omega_err[i]));
		fprintf(fp, " %.5e %.5e", imag(TP_tr->chi_omega[i]), imag(TP_tr->chi_omega_err[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	#endif // CHI_TR
	
}

void print_Gf0(int num)
{
	FILE *fp;
	char filename[128];
	
	complex<double> Gf0_omega[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		double omega_f = (double)(2*i+1) * M_PI / prm.beta;
		Gf0_omega[i] = 1./ ( IMAG * omega_f - prm.ef[0] - Delta_omega[i] );
	}
	
	double Gf0_tau[N_TAU+1];
	fft_fermion_radix2_omega2tau(Gf0_tau, Gf0_omega, prm.beta, N_TAU, 1);
	Gf0_tau[N_TAU] = -1.0 - Gf0_tau[0];
	
	sprintf(filename, DATA_DIR "%02d-Gf0_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		double tau = (double)i * prm.beta / (double)N_TAU;
		double chi0 = Gf0_tau[i] * Gf0_tau[N_TAU-i];
		fprintf(fp, "%.5e", tau);
		fprintf(fp, " %.6e", Gf0_tau[i]);
		if( i <= N_TAU/2 ){
			fprintf(fp, " %.6e", chi0);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	complex<double> G0_omega[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++)  G0_omega[i] = Delta_omega[i] / prm_V_sqr;
	double G0_tau[N_TAU+1];
	fft_fermion_radix2_omega2tau(G0_tau, G0_omega, prm.beta, N_TAU, 1);
	G0_tau[N_TAU] = -1.0 - G0_tau[0];
	
	sprintf(filename, DATA_DIR "%02d-G0_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		double tau = (double)i * prm.beta / (double)N_TAU;
		fprintf(fp, "%.5e", tau);
// 		fprintf(fp, " %.6e", G0_tau[i]);
		fprintf(fp, " %.6e", G0_tau[i]*prm_V_sqr);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
}


void print_D0(int num)
{
	FILE *fp;
	char filename[128];
	
	double D0_tau[N_TAU+1];
	fft_boson_radix2_omega2tau(D0_tau, D0_omega, prm.beta, N_TAU);
	D0_tau[N_TAU] = D0_tau[0];
	
	sprintf(filename, DATA_DIR "%02d-D0_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		double tau = (double)i * prm.beta / (double)N_TAU;
// 		fprintf(fp, "%.6e %.6e", tau, D0_tau[i]);
		fprintf(fp, "%.6e %.6e", tau, -D0_tau[i] * prm_gsqr_ch);
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	complex<double> check_D0_omega[N_TAU/2+1];
	fft_boson_radix2_tau2omega(D0_tau, check_D0_omega, prm.beta, N_TAU);
	
	sprintf(filename, DATA_DIR "%02d-D0_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU/2; i++){
		fprintf(fp, "%d", i);
// 		fprintf(fp, " %.6e %.6e", real(D0_omega[i]), real(check_D0_omega[i]));
		fprintf(fp, " %.6e", real(D0_omega[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	for(int i=0; i<=N_TAU/2; i++){
		if( abs(D0_omega[i] - check_D0_omega[i]) > 1e-10 ){
			printf("*** D0 i=%d %.5e\n", i, abs(D0_omega[i] - check_D0_omega[i]));
			exit(0);
		}
	}
	
	
	if( prm_dos_boson == 0 ){
		double K_tau_exact[N_TAU+1];
	// 	K_tau_exact[0] = K_tau_exact[N_TAU] = 0;
		for(int i=0; i<=N_TAU; i++){
			double tau = (double)i * prm.beta / (double)N_TAU;
			double t1 = prm.beta / 2. - tau;
			double t2 = prm.beta / 2.;
			
			K_tau_exact[i] = -( cosh(t1 * prm_w0) - cosh(t2 * prm_w0) ) / sinh( t2 * prm_w0 );
// 			K_tau_exact[i] *= pow(prm_g / prm_w0, 2);
			K_tau_exact[i] *= prm_gsqr_z / pow(prm_w0, 2);
		}
		
		
		double K_tau[N_TAU+1];
		K_tau[0] = K_tau[N_TAU] = 0;
		for(int i=1; i<N_TAU; i++){
			double tau = (double)i * prm.beta / (double)N_TAU;
			double beta_tau = prm.beta - tau;
			double t1 = prm.beta/2. - tau;
			double t2 = prm.beta/2.;
			
			K_tau[i] = -0.5 * real(D0_omega[0]) * tau * beta_tau;
	// 		K_tau[i] = real(D0_omega[0]) * 0.5*prm.beta * (0.5*prm.beta - tau) / 4.;
	// 		K_tau[i] = 0;
			
	// 		double s = -1;
	// 		for(int j=1; j<=N_TAU/2; j++){
	// 			double omega_b = double(2*j) * M_PI / prm.beta;
	// 			K_tau[i] += -s * 2.* real(D0_omega[j]) / (omega_b*omega_b) * ( cos(t1 * omega_b) - cos(t2 * omega_b) );
	// 			s = -s;
	// 		}
			
	// 		for(int j=1; j<=N_TAU/2; j++){
	// 			double omega_b = double(2*j) * M_PI / prm.beta;
	// 			K_tau[i] += -real(D0_omega[j]) / (omega_b*omega_b) * ( cos( tau * omega_b) + cos( beta_tau * omega_b) -2. );
	// 		}
			
			for(int j=1; j<=N_TAU/2; j++){
				double omega_b = double(2*j) * M_PI / prm.beta;
				K_tau[i] += 2.* real(D0_omega[j]) / (omega_b*omega_b) * ( 1.- cos( tau * omega_b) );
			}
			
			K_tau[i] *= prm_gsqr_z / prm.beta;
		}
		
		for(int i=0; i<=N_TAU/2; i++){
			if( abs(K_tau[i] - K_tau_exact[i]) > 1e-8 ){
				printf("*** K0 i=%d %.5e\n", i, abs(K_tau[i] - K_tau_exact[i]));
				exit(0);
			}
		}
		
		sprintf(filename, DATA_DIR "%02d-K0_tau.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<=N_TAU; i++){
			double tau = (double)i * prm.beta / (double)N_TAU;
			fprintf(fp, "%.6e %.6e %.6e", tau, K_tau[i], K_tau_exact[i]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
}


void print_vertex(int num)
{
	FILE *fp;
	char filename[128];
	
	sprintf(filename, DATA_DIR "%02d-vx-Gf.dat", num);
	fp=fopen(filename, "w");
	for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
		fprintf(fp, "%d", (2*(iw-N_VX1)+1));
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(VX->G[s][iw]), imag(VX->G[s][iw]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx-self.dat", num);
	fp=fopen(filename, "w");
	for(int iw=0; iw<2*N_VX1+N_VX2; iw++){
		fprintf(fp, "%d", (2*(iw-N_VX1)+1));
		for(int s=0; s<2; s++){
			fprintf(fp, " %.6e %.6e", real(VX->self[s][iw]), imag(VX->self[s][iw]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	complex<double> suscep_sp[N_VX2];
	complex<double> suscep_ch[N_VX2];
	for(int nu=0; nu<N_VX2; nu++){
		complex<double> diag = VX->suscep_lo[0][0][nu] + VX->suscep_lo[1][1][nu];
		complex<double> offd = VX->suscep_lo[0][1][nu] + VX->suscep_lo[1][0][nu];
		suscep_ch[nu] = diag + offd;
		suscep_sp[nu] = diag - offd;
	}
	
	sprintf(filename, DATA_DIR "%02d-vx-suscep.dat", num);
	fp=fopen(filename, "w");
	for(int nu=0; nu<N_VX2; nu++){
		fprintf(fp, "%d", nu);
		fprintf(fp, " %.6e", real(suscep_sp[nu]));
		fprintf(fp, " %.6e", real(suscep_ch[nu]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	// =======================================================
	
	complex<double> (*vx_sp)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx_ch)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx_pm)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx2_sp)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx2_ch)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx2_pm)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			for(int nu=0; nu<N_VX2; nu++){
				complex<double> diag, offd;
				diag = (VX->gamma_lo[0][0][iw1][iw2][nu] + VX->gamma_lo[1][1][iw1][iw2][nu]) / 2.;
				offd = (VX->gamma_lo[0][1][iw1][iw2][nu] + VX->gamma_lo[1][0][iw1][iw2][nu]) / 2.;
				vx_ch[iw1][iw2][nu] = diag + offd;
				vx_sp[iw1][iw2][nu] = diag - offd;
				vx_pm[iw1][iw2][nu] = (VX->gamma_tr[0][iw1][iw2][nu] + VX->gamma_tr[1][iw1][iw2][nu]) / 2.;
// 				vx_uu[iw1][iw2][nu] = diag;
// 				vx_ud[iw1][iw2][nu] = offd;
				
				diag = (VX->gamma2_lo[0][0][iw1][iw2][nu] + VX->gamma2_lo[1][1][iw1][iw2][nu]) / 2.;
				offd = (VX->gamma2_lo[0][1][iw1][iw2][nu] + VX->gamma2_lo[1][0][iw1][iw2][nu]) / 2.;
				vx2_ch[iw1][iw2][nu] = diag + offd;
				vx2_sp[iw1][iw2][nu] = diag - offd;
				vx2_pm[iw1][iw2][nu] = (VX->gamma2_tr[0][iw1][iw2][nu] + VX->gamma2_tr[1][iw1][iw2][nu]) / 2.;
// 				vx2_uu[iw1][iw2][nu] = diag;
// 				vx2_ud[iw1][iw2][nu] = offd;
				
				if( SELF_VERTEX_ISO ){
					vx_sp[iw1][iw2][nu] = (1./3.) * vx_sp[iw1][iw2][nu] + (2./3.) * vx_pm[iw1][iw2][nu];
					vx2_sp[iw1][iw2][nu] = (1./3.) * vx2_sp[iw1][iw2][nu] + (2./3.) * vx2_pm[iw1][iw2][nu];
				}
			}
		}
	}
	
	// asymptotic form at high frequencies
	complex<double> (*vx_sp_asympt)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	complex<double> (*vx_ch_asympt)[2*N_VX1][N_VX2] = new complex<double> [2*N_VX1][2*N_VX1][N_VX2];
	{
		double chi_uu[N_TP2+1], chi_ud[N_TP2+1], chi_pm[N_TP2+1];
		for(int i=0; i<=N_TP2; i++){
			chi_uu[i] = real(TP->chi_ch_omega[i] + TP->chi_sp_omega[i]) / 4.;
			chi_ud[i] = real(TP->chi_ch_omega[i] - TP->chi_sp_omega[i]) / 4.;
			chi_pm[i] = real(TP->chi_sp_omega[i]) / 2.;
		}
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			for(int iw2=0; iw2<2*N_VX1; iw2++){
				for(int nu=0; nu<N_VX2; nu++){
					complex<double> asympt_uu = prm.U*prm.U * (chi_uu[abs(iw1-iw2)] - chi_uu[nu]);
					complex<double> asympt_ud = prm.U*prm.U * (chi_pm[abs(iw1-iw2)] - chi_ud[nu]) + prm.U;
					vx_ch_asympt[iw1][iw2][nu] = asympt_uu + asympt_ud;
					vx_sp_asympt[iw1][iw2][nu] = asympt_uu - asympt_ud;
	// 				vx_uu_asympt[iw1][iw2][nu] = asympt_uu;
	// 				vx_ud_asympt[iw1][iw2][nu] = asympt_ud;
				}
			}
		}
	}
	
	
	
	// pointers to arrays of size [2*N_VX1][2*N_VX1][N_VX2]
	int n_vx = 5;
	complex<double> (*vx[n_vx])[2*N_VX1][N_VX2];
	vx[0] = vx_sp;
	vx[1] = vx_ch;
// 	vx[0] = vx_uu;
// 	vx[1] = vx_ud;
	vx[2] = vx_pm;
	vx[3] = vx_sp_asympt;
	vx[4] = vx_ch_asympt;
// 	vx[3] = vx_uu_asympt;
// 	vx[4] = vx_ud_asympt;
	
	// -------------------------------------------------------
	// for contour plot
	
	sprintf(filename, DATA_DIR "%02d-vx4-gamma0.dat", num);
	fp=fopen(filename, "w");
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			int nu=0;
			fprintf(fp, "%d %d", (2*(iw1-N_VX1)+1), (2*(iw2-N_VX1)+1));
			for(int i=0; i<n_vx; i++){
				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), imag(vx[i][iw1][iw2][nu]));
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx4-gamma1.dat", num);
	fp=fopen(filename, "w");
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			int nu=1;
			fprintf(fp, "%d %d", (2*(iw1-N_VX1)+1), (2*(iw2-N_VX1)+1));
			for(int i=0; i<n_vx; i++){
				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), imag(vx[i][iw1][iw2][nu]));
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx4-gamma2.dat", num);
	fp=fopen(filename, "w");
	for(int iw1=0; iw1<2*N_VX1; iw1++){
		for(int iw2=0; iw2<2*N_VX1; iw2++){
			int nu=N_VX2-1;
			fprintf(fp, "%d %d", (2*(iw1-N_VX1)+1), (2*(iw2-N_VX1)+1));
			for(int i=0; i<n_vx; i++){
				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), imag(vx[i][iw1][iw2][nu]));
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	// -------------------------------------------------------
	// nu-dependence
	
// 	int wf[][2] = {
// 		{N_VX1, N_VX1},
// 		{N_VX1*2-1, N_VX1*2-1},
// 		{N_VX1*2-1, N_VX1},
// 		{N_VX1*2-1, 0},
// 	};
// 	sprintf(filename, DATA_DIR "%02d-vx4-nu.dat", num);
// 	fp=fopen(filename, "w");
// 	fprintf(fp, "# (w1,w2) =");
// 	for(int j=0; j<4; j++)  fprintf(fp, " (%d,%d)", 2*(wf[j][0]-N_VX1)+1, 2*(wf[j][1]-N_VX1)+1);
// 	fprintf(fp, "\n");
// 	for(int nu=N_VX2-1; nu>0; nu--){  // nu<0
// 		fprintf(fp, "%d", -nu);
// 		for(int j=0; j<4; j++){
// 			int iw1 = 2*N_VX1 - wf[j][0] - 1;
// 			int iw2 = 2*N_VX1 - wf[j][1] - 1;
// 			for(int i=0; i<n_vx; i++){
// 				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), -imag(vx[i][iw1][iw2][nu]));
// 			}
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	for(int nu=0; nu<N_VX2; nu++){  // nu>0
// 		fprintf(fp, "%d", nu);
// 		for(int j=0; j<4; j++){
// 			int iw1 = wf[j][0];
// 			int iw2 = wf[j][1];
// 			for(int i=0; i<n_vx; i++){
// 				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), imag(vx[i][iw1][iw2][nu]));
// 			}
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	
	// -------------------------------------------------------
	// w1 scan
	
	int wb[4] = {0, 1, 2, N_VX2-1};
// 	int wb[4] = {0, 2, 10, N_VX2-1};  // Fig.10
// 	int wb[4] = {0, 5, 10, N_VX2-1};  // Fig.11
	{
		int iw2 = N_VX1+N_VX1/2;
// 		int iw2 = N_VX1+4;  // Fig.10
// 		int iw2 = N_VX1-1;  // Fig.11
		sprintf(filename, DATA_DIR "%02d-vx4-gamma-w.dat", num);
		fp=fopen(filename, "w");
		fprintf(fp, "# w2 = %d\n", (2*(iw2-N_VX1)+1));
		fprintf(fp, "# nu =");
		for(int j=0; j<4; j++)  fprintf(fp, " %d", wb[j]);
		fprintf(fp, "\n");
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			fprintf(fp, "%d", (2*(iw1-N_VX1)+1));
			for(int j=0; j<4; j++){
				int nu = wb[j];
				for(int i=0; i<n_vx; i++){
					fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw2][nu]), imag(vx[i][iw1][iw2][nu]));
				}
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
		
		
		sprintf(filename, DATA_DIR "%02d-vx4-gamma-w-test.dat", num);
		fp=fopen(filename, "w");
		for(int iw1=0; iw1<2*N_VX1; iw1++){
			fprintf(fp, "%d", (2*(iw1-N_VX1)+1));
			for(int j=0; j<4; j++){
				int nu = wb[j];
				fprintf(fp, " %.6e %.6e", real(vx2_sp[iw1][iw2][nu]), imag(vx2_sp[iw1][iw2][nu]));
				fprintf(fp, " %.6e %.6e", real(vx2_ch[iw1][iw2][nu]), imag(vx2_ch[iw1][iw2][nu]));
				fprintf(fp, " %.6e %.6e", real(vx2_pm[iw1][iw2][nu]), imag(vx2_pm[iw1][iw2][nu]));
// 				fprintf(fp, " %.6e %.6e", real(vx2_uu[iw1][iw2][nu]), imag(vx2_uu[iw1][iw2][nu]));
// 				fprintf(fp, " %.6e %.6e", real(vx2_ud[iw1][iw2][nu]), imag(vx2_ud[iw1][iw2][nu]));
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
	
	// -------------------------------------------------------
	// diagonal w1==w2
	
// 	sprintf(filename, DATA_DIR "%02d-vx4-diag.dat", num);
// 	fp=fopen(filename, "w");
// 	fprintf(fp, "# nu =");
// 	for(int j=0; j<4; j++)  fprintf(fp, " %d", wb[j]);
// 	fprintf(fp, "\n");
// 	for(int iw1=0; iw1<2*N_VX1; iw1++){
// 		fprintf(fp, "%d", (2*(iw1-N_VX1)+1));
// 		for(int j=0; j<4; j++){
// 			int nu = wb[j];
// 			for(int i=0; i<n_vx; i++){
// 				fprintf(fp, " %.6e %.6e", real(vx[i][iw1][iw1][nu]), imag(vx[i][iw1][iw1][nu]));
// 			}
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	
	// =======================================================
	
	complex<double> (*Gthree_sp)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	complex<double> (*Gthree_ch)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	complex<double> (*lambda_sp)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	complex<double> (*lambda_ch)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	complex<double> (*lambda2_sp)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	complex<double> (*lambda2_ch)[N_VX2] = new complex<double> [2*N_VX1][N_VX2];
	for(int iw=0; iw<2*N_VX1; iw++){
		for(int nu=0; nu<N_VX2; nu++){
			complex<double> diag, offd;
			diag = (VX->Gthree_lo[0][0][iw][nu] + VX->Gthree_lo[1][1][iw][nu]) / 2.;
			offd = (VX->Gthree_lo[0][1][iw][nu] + VX->Gthree_lo[1][0][iw][nu]) / 2.;
			Gthree_ch[iw][nu] = diag + offd;
			Gthree_sp[iw][nu] = diag - offd;
// 					Gthree_ch[iw][nu] = diag;
// 					Gthree_sp[iw][nu] = offd;
			
			diag = (VX->lambda_lo[0][0][iw][nu] + VX->lambda_lo[1][1][iw][nu]) / 2.;
			offd = (VX->lambda_lo[0][1][iw][nu] + VX->lambda_lo[1][0][iw][nu]) / 2.;
			lambda_ch[iw][nu] = diag + offd;
			lambda_sp[iw][nu] = diag - offd;
// 					lambda_ch[iw][nu] = diag;
// 					lambda_sp[iw][nu] = offd;
// 					lambda_ch[iw][nu] /= suscep_ch[nu] / 2.;
// 					lambda_sp[iw][nu] /= suscep_sp[nu] / 2.;
			
			diag = (VX->lambda_lo_test[0][0][iw][nu] + VX->lambda_lo_test[1][1][iw][nu]) / 2.;
			offd = (VX->lambda_lo_test[0][1][iw][nu] + VX->lambda_lo_test[1][0][iw][nu]) / 2.;
			lambda2_ch[iw][nu] = diag + offd;
			lambda2_sp[iw][nu] = diag - offd;
// 					lambda2_ch[iw][nu] = diag;
// 					lambda2_sp[iw][nu] = offd;
// 					lambda2_ch[iw][nu] /= suscep_ch[nu] / 2.;
// 					lambda2_sp[iw][nu] /= suscep_sp[nu] / 2.;
		}
	}
	
	int nu_plot[4] = {0, 1, N_VX2/2, N_VX2-1};
// 	int nu_plot[4] = {0, 1, 4, 9};
// 	sprintf(filename, DATA_DIR "%02d-vx3-chi.dat", num);
// 	fp=fopen(filename, "w");
// 	fprintf(fp, "# nu =");
// 	for(int j=0; j<4; j++)  fprintf(fp, " %d", nu_plot[j]);
// 	fprintf(fp, "\n");
// 	for(int nu=0; nu<N_VX2; nu++){
// 		for(int iw=0; iw<2*N_VX1; iw++){
// 			fprintf(fp, "%d %d", nu, (2*(iw-N_VX1)+1));
// 			fprintf(fp, " %.6e %.6e", real(Gthree_sp[iw][nu]), imag(Gthree_sp[iw][nu]));
// 			fprintf(fp, " %.6e %.6e", real(Gthree_ch[iw][nu]), imag(Gthree_ch[iw][nu]));
// 			fprintf(fp, "\n");
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	
// 	sprintf(filename, DATA_DIR "%02d-vx3-chi-w.dat", num);
// 	fp=fopen(filename, "w");
// 	fprintf(fp, "# nu =");
// 	for(int j=0; j<4; j++)  fprintf(fp, " %d", nu_plot[j]);
// 	fprintf(fp, "\n");
// 	for(int iw=0; iw<2*N_VX1; iw++){
// 		fprintf(fp, "%d", (2*(iw-N_VX1)+1));
// 		for(int j=0; j<4; j++){
// 			int nu = nu_plot[j];
// 			fprintf(fp, " %.6e %.6e", real(Gthree_sp[iw][nu]), imag(Gthree_sp[iw][nu]));
// 			fprintf(fp, " %.6e %.6e", real(Gthree_ch[iw][nu]), imag(Gthree_ch[iw][nu]));
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx3-lambda.dat", num);
	fp=fopen(filename, "w");
	for(int nu=0; nu<N_VX2; nu++){
		for(int iw=0; iw<2*N_VX1; iw++){
			fprintf(fp, "%d %d", nu, (2*(iw-N_VX1)+1));
			fprintf(fp, " %.6e %.6e", real(lambda_sp[iw][nu]), imag(lambda_sp[iw][nu]));
			fprintf(fp, " %.6e %.6e", real(lambda_ch[iw][nu]), imag(lambda_ch[iw][nu]));
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx3-lambda-w.dat", num);
	fp=fopen(filename, "w");
	fprintf(fp, "# nu =");
	for(int j=0; j<4; j++)  fprintf(fp, " %d", nu_plot[j]);
	fprintf(fp, "\n");
	for(int iw=0; iw<2*N_VX1; iw++){
		fprintf(fp, "%d", (2*(iw-N_VX1)+1));
		for(int j=0; j<4; j++){
			int nu = nu_plot[j];
			fprintf(fp, " %.6e %.6e", real(lambda_sp[iw][nu]), imag(lambda_sp[iw][nu]));
			fprintf(fp, " %.6e %.6e", real(lambda_ch[iw][nu]), imag(lambda_ch[iw][nu]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-vx3-lambda-w-test.dat", num);
	fp=fopen(filename, "w");
	for(int iw=0; iw<2*N_VX1; iw++){
		fprintf(fp, "%d", (2*(iw-N_VX1)+1));
		for(int j=0; j<4; j++){
			int nu = nu_plot[j];
			fprintf(fp, " %.6e %.6e", real(lambda2_sp[iw][nu]), imag(lambda2_sp[iw][nu]));
			fprintf(fp, " %.6e %.6e", real(lambda2_ch[iw][nu]), imag(lambda2_ch[iw][nu]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
// 	sprintf(filename, DATA_DIR "%02d-vx3-lambda-nu.dat", num);
// 	fp=fopen(filename, "w");
// 	for(int nu=0; nu<N_VX2; nu++){
// 		fprintf(fp, "%d", nu);
// 		int w[2] = {N_VX1+N_VX1/2, N_VX1/2};
// 		for(int j=0; j<2; j++){
// 			int iw = w[j];
// 			fprintf(fp, " %.6e %.6e", real(lambda_sp[iw][nu]), imag(lambda_sp[iw][nu]));
// 			fprintf(fp, " %.6e %.6e", real(lambda_ch[iw][nu]), imag(lambda_ch[iw][nu]));
// 		}
// 		fprintf(fp, "\n");
// 	}
// 	fclose(fp);
// 	printf(" '%s'\n", filename);
// 	
	
	delete [] vx_sp;
	delete [] vx_ch;
	delete [] vx_pm;
	delete [] vx2_sp;
	delete [] vx2_ch;
	delete [] vx2_pm;
	delete [] vx_sp_asympt;
	delete [] vx_ch_asympt;
	
	delete [] Gthree_sp;
	delete [] Gthree_ch;
	delete [] lambda_sp;
	delete [] lambda_ch;
	delete [] lambda2_sp;
	delete [] lambda2_ch;
}


static double func_energy_kin(double x, void *params)
{
	double *param = (double *)params;
	double gamma = param[0];
	double y = param[1];
	
	return( pow(x, gamma) / pow( y*y + x*x, 2 ) * ( x*(x*x - y*y) ) );
}
static double func_energy_int(double x, void *params)
{
	double *param = (double *)params;
	double gamma = param[0];
	double y = param[1];
	
	return( pow(x, gamma) / ( y*y + x*x ) * ( -2.*x ) );
}

// static double func_energy_re(double x, void *params)
// {
// 	double *param = (double *)params;
// 	double gamma = param[0];
// 	double y = param[1];
// 	
// 	return( pow(x, gamma) / pow( y*y + x*x, 2 ) * ( -x*(x*x + 3.*y*y) ) );
// }
// static double func_energy_im(double x, void *params)
// {
// 	double *param = (double *)params;
// 	double gamma = param[0];
// 	double y = param[1];
// 	
// 	return( pow(x, gamma) / pow( y*y + x*x, 2 ) * ( -2.*y*y*y ) );
// }

static double eval_energy_sub(complex<double> *chi, double *fac, double g_sqr)
{
	double energy = fac[0] * real(chi[0]);
	for(int i=1; i<=N_TAU/2; i++){
		energy += 2.* fac[i] * real(chi[i]);
	}
	energy *= g_sqr / prm.beta;
	
	return energy;
}

// DOS \propto omega^{gamma}  (with cutoff)
static void eval_energy(double beta, double cutoff, double gamma)
{
// 	complex<double> fac[N_TAU/2+1];
	double fac_kin[N_TAU/2+1];
	double fac_int[N_TAU/2+1];
	{
	// 	D0_omega[0] = 1./ gamma;
		
		gsl_function F;
		double params[2] = {gamma};
		F.params = params;
		
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
		int key = 6;
		double result, error;
		
		for(int i=0; i<=N_TAU/2; i++){
			double omega_b = double(2*i) * M_PI / beta;
			params[1] = omega_b / cutoff;
			
// 			F.function = &func_energy_re;
// 			gsl_integration_qag (&F, 0, 1, 0, 1e-7, 1000, key, w, &result_re, &error);
			
// 			F.function = &func_energy_im;
// 			gsl_integration_qag (&F, 0, 1, 0, 1e-7, 1000, key, w, &result_im, &error);
			
// 			fac[i] = complex<double>(result_re, result_im);
// 			fac[i] = result_re;
			
			F.function = &func_energy_kin;
			gsl_integration_qag (&F, 0, 1, 0, 1e-7, 1000, key, w, &result, &error);
			fac_kin[i] = result * (gamma+1.) / cutoff;
			
			F.function = &func_energy_int;
			gsl_integration_qag (&F, 0, 1, 0, 1e-7, 1000, key, w, &result, &error);
			fac_int[i] = result * (gamma+1.) / cutoff;
		}
		gsl_integration_workspace_free (w);
	}
// 	printf("\n%lf %lf\n", real(fac[0]), imag(fac[0]));
// 	for(int i=0; i<=N_TAU/2; i++)  printf("%lf %lf\n", fac[i]);
// 	exit(0);
	
	
	if( prm_gsqr_z == prm_gsqr_xy ){
		complex<double> *chi = TP->chi_sp_omega;
		double stat_suscep_err = PQ->stat_suscep_sp_err;
		if( flag_tp == 0 ){
			chi = TP->chi_pm_omega;
			stat_suscep_err = PQ->stat_suscep_pm_err;
		}
		
		PQ_energy_kin = 3.* eval_energy_sub(chi, fac_kin, prm_gsqr_z);
		PQ_energy_int = 3.* eval_energy_sub(chi, fac_int, prm_gsqr_z);
		
// 		double energy_kin_err
// 		 = fac_kin[0] * PQ->stat_suscep_sp_err * 3.* prm_gsqr_z / prm.beta;
// 		double energy_int_err
// 		 = fac_int[0] * PQ->stat_suscep_sp_err * 3.* prm_gsqr_z / prm.beta;
// 		
// 		printf("\n Internal energy\n  %lf\n", PQ_energy_tot);
// 		printf("  %lf %lf\n", PQ_energy_kin, energy_kin_err);
// 		printf("  %lf %lf\n", PQ_energy_int, energy_int_err);
		
		PQ_energy_tot_err = (fac_kin[0] + fac_int[0]) * stat_suscep_err
		 * 3.* prm_gsqr_z / prm.beta;
	}
	else{
		if( flag_tp ){
			PQ_energy_kin = eval_energy_sub(TP->chi_sp_omega, fac_kin, prm_gsqr_z);
			PQ_energy_int = eval_energy_sub(TP->chi_sp_omega, fac_int, prm_gsqr_z);
			
			PQ_energy_kin += 2.* eval_energy_sub(TP->chi_pm_omega, fac_kin, prm_gsqr_xy);
			PQ_energy_int += 2.* eval_energy_sub(TP->chi_pm_omega, fac_int, prm_gsqr_xy);
			
			PQ_energy_tot_err = (fac_kin[0] + fac_int[0]) * PQ->stat_suscep_sp_err
			 * prm_gsqr_z / prm.beta;
			PQ_energy_tot_err += (fac_kin[0] + fac_int[0]) * PQ->stat_suscep_pm_err
			 * 2.* prm_gsqr_xy / prm.beta;
		}
		else{
			PQ_energy_kin = PQ_energy_int = 0;
		}
	}
	
	double energy_h = -prm_h * PQ->occup_dif / 2.;
	PQ_energy_tot = PQ_energy_kin + PQ_energy_int + energy_h;
}

