/*

Continuous-time quantum Monte Carlo
for the impurity Anderson model (hybridization expansion)

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

static const char tag[64] = "v1.12";

#include "hyb_qmc.h"
#include "green_func_0.h"
#include "dmft.h"
#include "fft.h"
#include "pade.h"


#define DMFT_ISO 1
// 0: inisotropic bath, 1: isotropic bath

#define DMFT_OUTPUT_ITER 0
// output cavity Green function at each iteraction

int DMFT_ITER=0;
double DMFT_UPDATE_RATE=1.0;

int N_PADE = MIN(N_TAU/2, 8192);
// Number of frequencies used for Pade approximation of Gf (N_PADE^2 memory)
// MAX: N_TAU/2 (all data),  MIN: 0 (do not evaluate)


struct hyb_qmc_params prm;
double prm_ave_n, prm_chem_pot;
double prm_V_sqr[N_S], prm_V_sqr_imp[N_S], prm_D, prm_f_disp=0;
double prm_E_f[N_S], prm_E_c[N_S];

struct num_mc n_mc;

struct phys_quant *PQ;
struct single_particle *SP;
struct two_particle **TP, *TP_sp, *TP_ch;
struct two_particle_tr *TP_tr;
double *TP_tau;
struct phonon_g *D;

complex<double> G0_omega[N_S][N_TAU/2];

// struct dmft_green_func G;
struct dmft_green_func *G;

double SP_c_number[N_S];

int N_T, i_tmesh;
double T_MIN, T_MAX;

int my_rank;


int i_program;
char str_program[3][20]={
	"T,H fixed",
	"T varies",
	"H varies"
};

int flag_tp;
char str_flag[2][100]={
	"single-particle Green function",  // 0
	"single-particle & two-particle Green function",  // 1
};

int flag_dmft;
char str_dmft[2][100]={
	"impurity problem",
	"DMFT",
};

int flag_ensemble;
char str_ensemble[2][100]={
	"grand canonical ensemble (chemical potential given)",
	"canonical ensemble (band filling given)"
};



// ============================================================================

void init_params()
{
	FILE *fp;
	if( (fp = fopen("hyb_qmc_param.init", "r")) == NULL ){
		printf("'hyb_qmc_param.init' not found\n");
		exit(0);
	}
	
	if(my_rank==0){
		printf("\nCT-QMC for the Anderson model (CT-HYB)  %s\n", tag);
		printf(" DOS = %d\n", DOS);
	}
	
	double data[16];
	
	if( read_data(fp, data) != 1 )  exit(0);
	i_program = (int)data[0];
	
	if( read_data(fp, data) != 1 )  exit(0);
	flag_tp = (int)data[0];
	
	switch( read_data(fp, data) ){
		case 4:
			prm_f_disp = data[3];
		case 3:
			DMFT_UPDATE_RATE = data[2];
		case 2:
			DMFT_ITER = (int)data[1];
		case 1:
			flag_dmft = (int)data[0];
			break;
		default:
			exit(0);
	}
	
	if( read_data(fp, data) != 2 )  exit(0);
	flag_ensemble = (int)data[0];
	if( flag_ensemble==0 ){  // mu given
		prm_chem_pot = data[1];
	}
	else{  // n given
		prm_ave_n = data[1];
		prm_chem_pot = 0;
	}
	
	double prm_epsilon_f, prm_H;
	
	double prm_U;
	if( fscanf(fp, "%lf %lf %lf %lf", &prm_D, &prm_V_sqr[0], &prm_epsilon_f, &prm_U) != 4 )  exit(0);
	if( fscanf(fp, "%lf %lf", &prm.beta, &prm_H) != 2 )  exit(0);
	if( isinf(prm_U) )  prm.UINF = 1;
	
	if(PHONON){
		if( fscanf(fp, "%lf %lf", &prm.w_ph, &prm.g) != 2 )  exit(0);
	}
	
	if(my_rank==0){
		printf("\n %d: %s\n", i_program, str_program[i_program]);
		printf(" %d: %s\n", flag_tp, str_flag[flag_tp]);
		printf(" %d: %s\n", flag_dmft, str_dmft[flag_dmft]);
		printf(" %d: %s\n", flag_ensemble, str_ensemble[flag_ensemble]);
		
		printf("\n D = %.4e\n", prm_D);
		printf(" V_sqr = %.4e\n", prm_V_sqr[0]);
		printf(" epsilon_f = %.4e\n", prm_epsilon_f);
		printf(" U = %.4e  (UINF = %d)\n", prm_U, prm.UINF);
		printf(" beta = %.4e\n", prm.beta);
		printf(" H = %.4e\n", prm_H);
		if(PHONON){
			printf(" w_ph = %.4e\n", prm.w_ph);
			printf(" g = %.4e\n", prm.g);
		}
		
		if( !flag_ensemble )  printf(" chem_pot = %.4e\n", prm_chem_pot);  // grand canonical
		else                  printf(" ave_nc = %.4e\n", prm_ave_n);  // canonical
		
		if( flag_dmft==1 ){
			printf(" DMFT_ITER = %d\n", DMFT_ITER);
			printf(" DMFT_UPDATE_RATE = %.2lf\n", DMFT_UPDATE_RATE);
			printf(" f_disp = %.3lf\n", prm_f_disp);
		}
		printf("\n");
	}
	
	for(int s=1; s<N_S; s++)  prm_V_sqr[s] = prm_V_sqr[0];
	for(int s=0; s<N_S; s++)  prm_V_sqr_imp[s] = prm_V_sqr[s];
	for(int s1=0; s1<N_S; s1++){
		for(int s2=0; s2<N_S; s2++)  prm.U[s1][s2] = prm_U;
	}
	
	double moment_f[N_S];
	for(int s=0; s<N_S; s++){
		moment_f[s] = -(double)(N_S-1) * 0.5 + (double)s;
// 		prm.ef[s] = prm_epsilon_f + prm_H * moment_f[s];
		prm_E_f[s] = prm_epsilon_f + prm_H * moment_f[s];
		
		prm.ef[s] = prm_E_f[s];
		
		prm_E_c[s] = 0;
		
		if(my_rank==0){
			printf(" ef[%d] = %.5lf   ec[%d] = %.5lf\n", s, prm_E_f[s], s, prm_E_c[s]);
		}
	}
	
	
	if( fscanf(fp, "%d %d %d %d", &n_mc.N_BIN, &n_mc.N_MSR, &n_mc.N_ADD, &n_mc.N_SHIFT) != 4 )  exit(0);
	if(my_rank==0){
		printf("\n N_BIN = %d  N_MSR = %d  (N_ADD = %d  N_SHIFT = %d)\n",
		 n_mc.N_BIN, n_mc.N_MSR, n_mc.N_ADD, n_mc.N_SHIFT);
	}
	
	if(i_program==1){
		if( fscanf(fp, "%d %d %lf %lf", &i_tmesh, &N_T, &T_MIN, &T_MAX) != 4 )  exit(0);
	}
	
	fclose(fp);
	
	
// 	#if PHONON
// 	for(int s1=0; s1<N_S; s1++){
// 		prm.ef[s1] += prm.g * prm.g / prm.w_ph;
// 		
// 		for(int s2=0; s2<N_S; s2++){
// 			prm.U[s1][s2] -= 2.0 * prm.g * prm.g / prm.w_ph;
// 		}
// 	}
// 	#endif // PHONON
	
}

// ============================================================================

void print_pq(int num)
{
	char filename[128];
	sprintf(filename, DATA_DIR "T_depend.dat");
	
	static FILE *fp=fopen(filename, "w");
	fp=fopen(filename, "a");
	
	fprintf(fp, "%.5e", 1./prm.beta);
	for(int s=0; s<N_S; s++)  fprintf(fp, " %.5e %.3e", SP[s].f_number, SP[s].f_number_err);
	fprintf(fp, " %.5e %.4e", PQ->stat_suscep_sp, PQ->stat_suscep_sp_err);
	fprintf(fp, " %.5e %.4e", PQ->stat_suscep_ch, PQ->stat_suscep_ch_err);
	// fprintf(fp, " %.5e %.4e", TP[0][1].chi_tau[0], TP[0][1].chi_tau_err[0]);
	for(int s=0; s<N_S; s++){
		fprintf(fp, " %.5e %.5e", real(SP[s].self_f[0]), imag(SP[s].self_f[0]));
		fprintf(fp, " %.5e %.5e", real(SP[s].self_f[1]), imag(SP[s].self_f[1]));
	}
	for(int s=0; s<N_S; s++)  fprintf(fp, " %.5e", SP_c_number[s]);
	if( flag_dmft ){
		fprintf(fp, " %.6e %.6e", prm_chem_pot, prm_ave_n);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.5e %.5e", G[s].n_f, G[s].n_c);
// 			fprintf(fp, " %.5e %.5e", real(G[s].F[0])*prm_V_sqr[s], imag(G[s].F[0])*prm_V_sqr[s]);
// 			fprintf(fp, " %.5e %.5e", real(G[s].F[1])*prm_V_sqr[s], imag(G[s].F[1])*prm_V_sqr[s]);
			fprintf(fp, " %.5e %.5e", real(G[s].self_c[0]), imag(G[s].self_c[0]));
			fprintf(fp, " %.5e %.5e", real(G[s].self_c[1]), imag(G[s].self_c[1]));
		}
	}
	if( PHONON ){
		fprintf(fp, " %.6e", D->occup);
	}
	
	fprintf(fp, "\n");
	fclose(fp);
	printf("\n '%s' updated\n", filename);
}

void print_statistics(int num)
{
	FILE *fp;
	char filename[40];
	
	sprintf(filename, DATA_DIR "%02d-stat.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_K * N_S; i++){
		fprintf(fp, "%d", i);
		if( PQ->Z_ktot[i] != 0. )  fprintf(fp, " %.5e", PQ->Z_ktot[i]);
		else  fprintf(fp, " ?");
		if(i < N_K){
			for(int s=0; s<N_S; s++){
				if( PQ->Z_k[s][i] != 0. )  fprintf(fp, " %.5e", PQ->Z_k[s][i]);
				else  fprintf(fp, " ?");
			}
		}
		fprintf(fp, "\n");
	}
	printf("\n '%s'\n", filename);
	fclose(fp);
}

void print_single_particle(int num)
{
	printf("\n f number                   (measured at tau=0)\n");
	for(int s=0; s<N_S; s++){
		printf("   %d: %lf +- %lf", s, SP[s].f_number, SP[s].f_number_err);
		printf("   (%lf +- %lf)\n", SP[s].f_number2, SP[s].f_number2_err);
	}
	printf(" tot: %lf +- %lf\n", PQ->occup_tot, PQ->occup_tot_err);
	
	FILE *fp;
	char filename[40];
	
	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}
	
	
	sprintf(filename, DATA_DIR "%02d-Gf_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.3e", tau[i]);
		
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", SP[s].Gf_tau[i], SP[s].Gf_tau_err[i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-Gf_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].Gf_omega[i]), imag(SP[s].Gf_omega[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	complex<double> i_omega_f[N_PADE];
	
	if(N_PADE){
		complex<double> Gf_pade[N_S][N_PADE];
		
		for(int i=0; i<N_PADE; i++){
			double omega_f = (double)(2*i+1) * M_PI / prm.beta;
			i_omega_f[i] = IMAG * omega_f;
		}
		
		for(int s=0; s<N_S; s++){
			pade_init(i_omega_f, SP[s].Gf_omega, Gf_pade[s], N_PADE);
		}
		
		sprintf(filename, DATA_DIR "%02d-Gf_pade.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<1000; i++){
			
			complex<double> w = -2 * prm_D + 4 * prm_D / 999. * (double)i;
			double delta0 = M_PI * prm_V_sqr_imp[0] / prm_D * 0.5;
			
			fprintf(fp, "%.6e", real(w));
			for(int s=0; s<N_S; s++){
				complex<double> temp_Gf = pade_calc(i_omega_f, Gf_pade[s], w, N_PADE);
				fprintf(fp, " %.6e %.6e", real(temp_Gf), imag(temp_Gf));
			}
			fprintf(fp, " %.6e",
			 -delta0 / (pow(real(w)-prm.ef[0]-prm.U[0][0]*0.5,2) + pow(delta0,2)) );
			
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
	
	
	sprintf(filename, DATA_DIR "%02d-self_f.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
// 	for(int i=0; i<100; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_f[i]), imag(SP[s].self_f[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	
	//
	// conduction electron Green function (hybridization is assumed to be local)
	//
	complex<double> Gc_omega[N_S][N_TAU/2];
	
	if(flag_dmft==0 || prm_f_disp==0){
		
		for(int s=0; s<N_S; s++){
			for(int i=0; i<N_TAU/2; i++){
				Gc_omega[s][i] = G0_omega[s][i] * ( 1.0 + G0_omega[s][i] * prm_V_sqr_imp[s] * SP[s].Gf_omega[i] );
			}
		}
		
		double Gc_tau[N_S][N_TAU];
		for(int s=0; s<N_S; s++){
			fft_fermion_radix2_omega2tau(Gc_tau[s], Gc_omega[s], prm.beta, N_TAU, 1);
			SP_c_number[s] = 1.0 + Gc_tau[s][0];
		}
		
		
		printf("\n occupation number  (total  c  f)\n");
		double tot_n=0, n_c=0, n_f=0;
		for(int s=0; s<N_S; s++){
			tot_n += (SP[s].f_number + SP_c_number[s]) / (double)N_S;
			n_c += SP_c_number[s] / (double)N_S;
			n_f += SP[s].f_number / (double)N_S;
			
			printf("  %3d: %lf  %lf  %lf\n", s, SP[s].f_number + SP_c_number[s], SP_c_number[s], SP[s].f_number);
		}
		printf("  ave: %lf  %lf  %lf\n", tot_n, n_c, n_f);
		
		
		sprintf(filename, DATA_DIR "%02d-Gc_omega.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<N_TAU/2; i++){
			fprintf(fp, "%d", i);
			for(int s=0; s<N_S; s++){
				fprintf(fp, " %.6e %.6e", real(Gc_omega[s][i]), imag(Gc_omega[s][i]));
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("\n '%s'\n", filename);
	}
	
	if(N_PADE){
		complex<double> Gc_pade[N_S][N_PADE];
		
		for(int s=0; s<N_S; s++){
			if(flag_dmft==0 || prm_f_disp==0){
				pade_init(i_omega_f, Gc_omega[s], Gc_pade[s], N_PADE);
			}
			else{
				pade_init(i_omega_f, G[s].Gc_perio, Gc_pade[s], N_PADE);
			}
		}
		
		sprintf(filename, DATA_DIR "%02d-Gc_pade.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<1000; i++){
			
			complex<double> w = -2 * prm_D + 4 * prm_D / 999. * (double)i;
			
			fprintf(fp, "%.6e", real(w));
			for(int s=0; s<N_S; s++){
				complex<double> temp_Gc = pade_calc(i_omega_f, Gc_pade[s], w, N_PADE);
				fprintf(fp, " %.6e %.6e", real(temp_Gc), imag(temp_Gc));
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
	
}


void print_two_particle(int num)
{
	printf("\n static susceptibility\n");
	printf("  spin   : %.5e +- %.4e\n", PQ->stat_suscep_sp, PQ->stat_suscep_sp_err);
	printf("  charge : %.5e +- %.4e\n", PQ->stat_suscep_ch, PQ->stat_suscep_ch_err);
	
// 	printf("\n static susceptibility (test)\n");
// 	printf("  spin   : %.5e\n", real(TP_sp->chi_omega[0]));
// 	printf("  charge : %.5e\n", real(TP_ch->chi_omega[0]));
	
	
	FILE *fp;
	char filename[40];
	
// 	double tau[N_TAU+1];
// 	for(int i=0; i<=N_TAU; i++){
// 		tau[i] = (double)i * prm.beta / (double)N_TAU;
// 	}
	
	//
	// non-interacting
	//
	double Gf0_tau[N_TAU+1];
	Gf_calc_fft(Gf0_tau, prm.beta, prm_D, prm.ef[0], prm.U[0][0], prm_V_sqr_imp[0]);
	Gf0_tau[N_TAU] = -1.0 - Gf0_tau[0];
	
	double chi0_tau[N_TAU+1];
	for(int i=0; i<=N_TAU/2; i++)  chi0_tau[i] = Gf0_tau[i] * Gf0_tau[N_TAU-i];
	
	sprintf(filename, DATA_DIR "%02d-chi0.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU/2; i++){
		double tau = (double)i * prm.beta / (double)N_TAU;
		fprintf(fp, "%.5e %.5e\n", tau, chi0_tau[i]);
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-chi_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TP; i++){
// 		double tau = (double)i * prm.beta / (double)(2*N_TP);
		fprintf(fp, "%.5e", TP_tau[i]);
		fprintf(fp, " %.5e %.4e", TP_sp->chi_tau[i], TP_sp->chi_tau_err[i]);
		fprintf(fp, " %.5e %.4e", TP_ch->chi_tau[i], TP_ch->chi_tau_err[i]);
		for(int s1=0; s1<N_S; s1++){
			for(int s2=0; s2<N_S; s2++){
				fprintf(fp, " %.5e %.4e", TP[s1][s2].chi_tau[i], TP[s1][s2].chi_tau_err[i]);
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	
	sprintf(filename, DATA_DIR "%02d-chi_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TP2; i++){
		fprintf(fp, "%d %.5e %.5e\n", i, real(TP_sp->chi_omega[i]), real(TP_ch->chi_omega[i]));
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
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
	
	
	//
	// phonon Green function
	//
	#if PHONON
	
	printf("\n phonon number: %.5e\n", D->occup);
	
	sprintf(filename, DATA_DIR "%02d-phonon_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TP2; i++){
		fprintf(fp, "%d", i);
		fprintf(fp, " %.5e", real(D->d_omega[i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	complex<double> d_pade[N_TP2];
	
	pade_init(i_omega_b, D->d_omega, d_pade, N_TP2);
	
	sprintf(filename, DATA_DIR "%02d-phonon_pade.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<1000; i++){
		complex<double> w = 1.0 * prm_D / 999. * (double)i;
		
		complex<double> temp_d = pade_calc(i_omega_b, d_pade, w, N_TP2);
		
		fprintf(fp, "%.6e %.6e %.6e\n", real(w),
		 real(temp_d), imag(temp_d));
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	#endif // PHONON
}


void print_G0_imp(int num)
{
	FILE *fp;
	char filename[40];
	
	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}
	
	
	double G0_tau[N_TAU+1];
	fft_fermion_radix2_omega2tau(G0_tau, G0_omega[0], prm.beta, N_TAU);
	
	sprintf(filename, DATA_DIR "%02d-Gc0_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU; i++){
		fprintf(fp, "%.4e %.6e\n", tau[i], G0_tau[i]);
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	
	sprintf(filename, DATA_DIR "%02d-Gc0_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
// 		double exact = -prm_V_sqr * atan(prm_D / omega_f[i]) / prm_D;
// 		double exact = -atan(prm_D / omega_f[i]) / prm_D;
		
// 		fprintf(fp, "%d %.6e %.6e\n", i, imag(G0_omega[0][i]), exact);
		fprintf(fp, "%d %.6e %.6e\n", i, real(G0_omega[0][i]), imag(G0_omega[0][i]));
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	if(N_PADE){
		complex<double> G0_pade[N_PADE], i_omega_f[N_PADE];
		
		for(int i=0; i<N_PADE; i++){
			double omega_f = (double)(2*i+1) * M_PI / prm.beta;
			i_omega_f[i] = IMAG * omega_f;
		}
		
		pade_init(i_omega_f, G0_omega[0], G0_pade, N_PADE);
		
		sprintf(filename, DATA_DIR "%02d-Gc0_pade.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<1000; i++){
			complex<double> temp = -1.1 * prm_D + 2.2 * prm_D / 999. * (double)i;
			complex<double> temp_G0 = pade_calc(i_omega_f, G0_pade, temp, N_PADE);
			
			fprintf(fp, "%.6e %.6e %.6e\n", real(temp), real(temp_G0), imag(temp_G0));
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
	
	//
	// Gf0
	//
	complex<double> Gf0_omega[N_TAU/2];
	
	sprintf(filename, DATA_DIR "%02d-Gf0_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		complex<double> i_omega_n = IMAG * (double)(2*i+1) * M_PI / prm.beta;
		Gf0_omega[i] = 1.0 / (i_omega_n - prm.ef[0] + prm_chem_pot - 0.5*prm.U[0][0]
		 - prm_V_sqr[0] * G0_omega[0][i]);
		
		fprintf(fp, "%d %.6e %.6e\n", i, real(Gf0_omega[i]), imag(Gf0_omega[i]));
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	double Gf0_tau[N_TAU+1];
	fft_fermion_radix2_omega2tau(Gf0_tau, Gf0_omega, prm.beta, N_TAU);
	
	sprintf(filename, DATA_DIR "%02d-Gf0_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU; i++){
		fprintf(fp, "%.4e %.6e\n", tau[i], Gf0_tau[i]);
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
/*	
	//
	// Gf0
	// for DOS==0 && symmetric condition
	//
	if(DOS==0){
		
		double Gf0_tau[N_S][N_TAU];
		for(int s=0; s<N_S; s++){
			Gf_calc_fft(Gf0_tau[s], prm.beta, prm_D, prm.ef[s], prm.U[0][0], prm_V_sqr_imp[0]);
		}
		
		sprintf(filename, DATA_DIR "%02d-Gf0_tau.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<N_TAU; i++){
			fprintf(fp, "%.4e", tau[i]);
			for(int s=0; s<N_S; s++){
				fprintf(fp, " %.6e", Gf0_tau[s][i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("\n '%s'\n", filename);
		
		
		complex<double> Gf0_omega[N_S][N_TAU/2];
		for(int s=0; s<N_S; s++){
			fft_fermion_radix2_tau2omega(Gf0_tau[s], Gf0_omega[s], prm.beta, N_TAU, 1);
		}
		
		sprintf(filename, DATA_DIR "%02d-Gf0_omega.dat", num);
		fp=fopen(filename, "w");
		for(int i=0; i<N_TAU/2; i++){
			fprintf(fp, "%d", i);
			for(int s=0; s<N_S; s++){
				fprintf(fp, " %.6e %.6e", real(Gf0_omega[s][i]), imag(Gf0_omega[s][i]));
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" '%s'\n", filename);
	}
*/	
	
}

void print_G0_lattice(int num)
{
	FILE *fp;
	char filename[40];
	
	int N = N_S;
	if(DMFT_ISO)  N = 1;
	
	
	sprintf(filename, DATA_DIR "%02d-Gc0_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<N; s++){
			fprintf(fp, " %.6e %.6e\n", real(G[s].Gc_perio[i]), imag(G[s].Gc_perio[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
	
	sprintf(filename, DATA_DIR "%02d-Gf0_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<N; s++){
			fprintf(fp, " %.6e %.6e\n", real(G[s].Gf_perio[i]), imag(G[s].Gf_perio[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);
	
	
	for(int s=0; s<N; s++){
		printf(" %d  n0=%.5lf  n0_f=%.5lf  n0_c=%.5lf\n", s, G[s].n_f + G[s].n_c, G[s].n_f, G[s].n_c);
	}
}


// ============================================================================

void dmft_iter(int num)
{
	int ns = DMFT_ISO ? 1 : N_S;
	
	for(int s=0; s<N_S; s++){
		prm_V_sqr_imp[s] = prm_V_sqr[s] + dmft_ads_Veff(prm_f_disp, prm_D);
	}
	if(my_rank==0){
		printf("\n---\n");
		printf("\n V_sqr_imp = %.4e\n", prm_V_sqr_imp[0]);
	}
	
	//
	// G0
	//
	complex<double> *p_self_f[N_S];
	
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU/2; i++)  SP[s].self_f[i] = prm.U[s][s] / N_S;  // Hartree term
		
		p_self_f[s] = SP[s].self_f;
	}
	
	dmft_ads1(G, p_self_f, ns, prm.beta, prm_D, prm_V_sqr, prm_E_f, prm_f_disp,
	 flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
	
	if(my_rank==0){
		print_G0_lattice(num);
	}
	
// 	void read_dmft_cav(int);
// 	read_dmft_cav(num);
	double mu_temp = 0;
	if( read_dmft_mu(G, ns, num, mu_temp) ){
		if( flag_ensemble )  prm_chem_pot = mu_temp;
		if(my_rank==0){
			printf("\nBath imported from an external file\n");
			printf(" Chemical potential: %.6e\n", prm_chem_pot);
		}
	}
	if( DMFT_ISO ){
		for(int s=1; s<N_S; s++)  G[s] = G[0];
	}
	
	
	for(int s=0; s<N_S; s++){
		for(int i=0; i<N_TAU/2; i++){
// 			G[s].Gc_cav_omega[i] = G0_omega[s][i];
			G0_omega[s][i] = G[s].Gc_cav_omega[i];
		}
	}
	
	
	FILE *fp_iter;
	char filename_iter[100];
	if(my_rank==0){
		sprintf(filename_iter, DATA_DIR "%02d-dmft-iter.dat", num);
		fp_iter=fopen(filename_iter, "w");
		fclose(fp_iter);
		printf("\n '%s'\n", filename_iter);
	}
	
// 	complex<double> pre_cav, pre_gc, pre_self, pre_tmat;
	
	for(int iter=0; iter<DMFT_ITER; iter++){
		
		char str_iter[100];
		if(my_rank==0){
			sprintf(str_iter, "\n---   iter=%d\n", iter);
			printf("%s", str_iter);
			hybqmc_fprint_log(str_iter);
		}
		
		double pre_chem_pot = prm_chem_pot;
		complex<double> pre_cav = G[0].Gc_cav_omega[0];
		complex<double> pre_gc = G[0].Gc_perio[0];
		complex<double> pre_F = G[0].F[0];
		complex<double> pre_tmat = G[0].tmat_perio[0];
		
		//
		// measure
		//
		for(int s=0; s<N_S; s++)  prm.ef[s] = prm_E_f[s] - prm_chem_pot;
		hybqmc_set_params(prm);
		
		hybqmc_set_G0(G0_omega, prm_V_sqr_imp);
		
		hybqmc_eval( iter == DMFT_ITER - 1 ? flag_tp : 0);
		
		
		if(my_rank==0){
			
			#if DMFT_ISO  // isotropic bath
			
			complex<double> ave_Gf_omega[N_TAU/2] = {0};
// 			double ave_f_number = 0;
// 			complex<double> ave_self_f[N_TAU/2] = {0};
			for(int s=0; s<N_S; s++){
// 				ave_f_number += SP[s].f_number / (double)N_S;
				
				for(int i=0; i<N_TAU/2; i++){
					ave_Gf_omega[i] += SP[s].Gf_omega[i] / (double)N_S;
					
// 					ave_self_f[i] += SP[s].self_f[i] / (double)N_S;
// 					ave_self_f[i] = 0;
					
// 					if( imag(ave_self_f[i]) > 0 )  ave_self_f[i] = real(ave_self_f[i]);
				}
			}
			
// 			dmft_ads0(G, ave_Gf_omega, prm.beta, prm_D, prm_V_sqr[0]);
// 			dmft_ads0(G[0], ave_Gf_omega, prm.beta, prm_D, prm_V_sqr[0], ave_f_number, flag_ensemble, prm_ave_n, prm_chem_pot);
			
// // 			complex<double> *p_self_f[1];
// 			p_self_f[0] = ave_self_f;
// 			
// 			dmft_ads1(G, p_self_f, 1, prm.beta, prm_D, prm_V_sqr, prm_E_f, prm_f_disp,
// 			 flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
			
			complex<double> *p_Gf[1];
			p_Gf[0] = ave_Gf_omega;
			
			dmft_ads2(G, p_Gf, 1, prm.beta, prm_D, prm_V_sqr, prm_f_disp,
			 flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
			
			
			#if DMFT_OUTPUT_ITER
			print_dmft_num(G, 1, num, iter);
			#endif // DMFT_OUTPUT_ITER
			
			for(int s=1; s<N_S; s++)  G[s] = G[0];
			
			
			#else // DMFT_ISO
			
// 			complex<double> *p_self_f[N_S];
// 			for(int s=0; s<N_S; s++)  p_self_f[s] = SP[s].self_f;
			
// 			dmft_ads1(G, p_self_f, N_S, prm.beta, prm_D, prm_V_sqr, prm_E_f, prm_f_disp,
// 			 flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
			
			complex<double> *p_Gf[N_S];
			for(int s=0; s<N_S; s++)  p_Gf[s] = SP[s].Gf_omega;
			
			dmft_ads2(G, p_Gf, N_S, prm.beta, prm_D, prm_V_sqr, prm_f_disp,
			 flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
			
			#if DMFT_OUTPUT_ITER
			print_dmft_num(G, N_S, num, iter);
			#endif // DMFT_OUTPUT_ITER
			
			#endif // DMFT_ISO
			
			
			double dif_chem_pot = prm_chem_pot - pre_chem_pot;
			complex<double> dif_cav = G[0].Gc_cav_omega[0] - pre_cav;
			complex<double> dif_gc = G[0].Gc_perio[0] - pre_gc;
			complex<double> dif_F = G[0].F[0] - pre_F;
			complex<double> dif_tmat = G[0].tmat_perio[0] - pre_tmat;
			
			fp_iter=fopen(filename_iter, "a");
			fprintf(fp_iter, "%d", iter);
			fprintf(fp_iter, " %.5e", dif_chem_pot);
			fprintf(fp_iter, " %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
			 real(dif_cav), imag(dif_cav), real(dif_gc), imag(dif_gc),
			 real(dif_F), imag(dif_F), real(dif_tmat), imag(dif_tmat));
			fclose(fp_iter);
			
		}
		
		
		#if HYB_QMC_MPI
		for(int s=0; s<N_S; s++){
			MPI_Bcast(G[s].Gc_cav_omega, N_TAU, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		// sizeof(complex<double>) * (N_TAU/2) = sizeof(double) * N_TAU
		
// 		MPI_Bcast(G.Gc_cav_omega, N_TAU, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif // HYB_QMC_MPI
		
		
		//
		// update
		//
		prm_chem_pot = DMFT_UPDATE_RATE * prm_chem_pot + (1.- DMFT_UPDATE_RATE) * pre_chem_pot;
		
		for(int s=0; s<N_S; s++){
			for(int i=0; i<N_TAU/2; i++){
// 				G0_omega[s][i] = G[s].Gc_cav_omega[i];
				G0_omega[s][i] = DMFT_UPDATE_RATE * G[s].Gc_cav_omega[i] + (1.- DMFT_UPDATE_RATE) * G0_omega[s][i];
				G[s].Gc_cav_omega[i] = G0_omega[s][i];
			}
		}
	}
	
	if(my_rank==0){
		printf("\n Chemical potential: %.6e\n", prm_chem_pot);
		
		printf(" Filling (ave.):   total   f   c\n");
		printf("   %.6e  %.6e  %.6e\n", prm_ave_n, G[0].n_f, G[0].n_c);
		
		
		double temp_E_f[N_S];
		for(int s=0; s<N_S; s++){
			temp_E_f[s] = prm_E_f[s];
			
// 			#if PHONON
// 			temp_E_f[s] += prm.g * prm.g / prm.w_ph;
// 			#endif // PHONON
		}
		
		
		print_dmft_mu(G, ns, num, prm_chem_pot);
		
		#if DMFT_ISO
// 		dmft_nk_ads(G[0].self_c, prm.beta, prm_D, prm_V_sqr[0], prm_chem_pot, num);
		
// 		complex<double> *p_self_f[1] = {SP[0].self_f};
		p_self_f[0] = SP[0].self_f;
		dmft_nk_ads(p_self_f, 1, prm.beta, prm_D, prm_V_sqr, temp_E_f, prm_f_disp, prm_E_c, prm_chem_pot, num);
		
		#else // DMFT_ISO
		dmft_nk_ads(p_self_f, N_S, prm.beta, prm_D, prm_V_sqr, temp_E_f, prm_f_disp, prm_E_c, prm_chem_pot, num);
		
		#endif // DMFT_ISO
		
		
		char str_iter[100];
		sprintf(str_iter, "\n---   iteration finish\n");
		printf("%s", str_iter);
		hybqmc_fprint_log(str_iter);
	}
	
}


void read_dmft_cav(int num)
{
	FILE *fp;
	char filename[40];
	
	sprintf(filename, "%02d-dmft.dat", num);
	
	if( (fp=fopen(filename, "r")) != NULL ){
		
		char str[1000];
		
		if(my_rank==0){
// 			printf("\n---\n");
			printf("\n'%s' opened\n", filename);
		}
		
		int n=0;  // data number
		
		while( fgets(str, sizeof(str), fp) != NULL ){
			if(str[0] != 35){  // 35=="#"
				
				#if DMFT_ISO  // isotropic bath
				
				int i_temp;
				double cav_re, cav_im;
				
				sscanf(str, "%d %lf %lf", &i_temp, &cav_re, &cav_im);
				G0_omega[0][n] = cav_re + IMAG * cav_im;
				
				for(int a=1; a<N_S; a++)  G0_omega[a][n] = G0_omega[0][n];
				
				#else // DMFT_ISO
				
				char *p_str;
				p_str = strtok(str, " ");
				p_str = strtok(NULL, " ");
				
				for(int a=0; a<N_S; a++){
					if(p_str != NULL){
// 						sscanf(p_str, "%lf %lf", &cav_re, &cav_im);
// 						G0_omega[a][n] = cav_re + IMAG * cav_im;
						
						double temp[8];
						for(int i=0; i<8; i++){
							sscanf(p_str, "%lf", &temp[i]);
							p_str = strtok(NULL, " ");
						}
						
						G0_omega[a][n] = temp[0] + IMAG * temp[1];
					}
					else if(a != 0){
						G0_omega[a][n] = G0_omega[0][n];
					}
					else{
						if(my_rank==0){
							printf("\n*** error : p_str == NULL\n");
						}
						exit(0);
					}
					
// 					for(int i=0; i<8; i++)  p_str = strtok(NULL, " ");
				}
				
				#endif // DMFT_ISO
				
				n++;
			}
		}
		
		fclose(fp);
		
		// test
// 		for(int i=0; i<n; i++){
// 			printf("%d", i);
// 			for(int a=0; a<N_F; a++)  printf(" %.5e %.5e", real(G0_omega[a][i]), imag(G0_omega[a][i]));
// 			printf("\n");
// 		}
		
		
		if(n != N_TAU/2){
			if(my_rank==0){
				printf("\n*** error : inconsistency in data number\n");
				printf("  n = %d,  N_TAU/2 = %d\n", n, N_TAU/2);
			}
			
			exit(0);
		}
		
		if(my_rank==0)  printf(" G0_omega[] has been set\n");
		
	}
	else{
		if(my_rank==0){
// 			printf("\n---\n");
			printf("\n'%s' does not exist\n", filename);
		}
	}
	
}


// ============================================================================


void main_sub(int num)
{
	clock_t time_start = clock();
	
// 	hybqmc_set_params(prm);
	
	if( !flag_dmft ){  // impurity problem
		for(int s=0; s<N_S; s++)  G0_omega_calc(G0_omega[s], prm.beta, prm_D);
		//	G0_omega_calc(G0_omega, N_S, prm.beta, prm_D, flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
		
		if(my_rank==0){
			print_G0_imp(num);
		}
		
		for(int s=0; s<N_S; s++)  prm.ef[s] = prm_E_f[s] - prm_chem_pot;
		hybqmc_set_params(prm);
		
		hybqmc_set_G0(G0_omega, prm_V_sqr_imp);
		
		hybqmc_eval(flag_tp);
	}
	else{  // DMFT
		dmft_iter(num);
	}
	
	
	if(my_rank==0){
		print_pq(num);
		print_statistics(num);
		print_single_particle(num);
		
		if(flag_tp){
			print_two_particle(num);
		}
		
		clock_t time_end = clock();
		print_time(time_start, time_end);
	// 	print_time(time_start, time_end, LOG_FILE);
	}
}

void main_T()
{
	clock_t time_start = clock();
	
	double temperature[N_T], T_energy[N_T], T_energy_err[N_T];
	
	mesh_init(i_tmesh, temperature, T_MIN, T_MAX, N_T, my_rank);
	
	
	for(int it=0; it<N_T; it++){
		
		if(my_rank==0){
			char str_t[100];
			sprintf(str_t, "\n===   it=%d  T=%lf\n", it, temperature[it]);
			printf("%s", str_t);
			hybqmc_fprint_log(str_t);
		}
		
		prm.beta = 1.0 / temperature[it];
		
		main_sub(it);
		
// 		T_energy[it] = PQ.energy;
// 		T_energy_err[it] = PQ.energy_err;
		
	}
	
	
	if(my_rank==0){
		printf("\n===\n");
		
	// 	sprintf(filename, DATA_DIR "T_thermo.dat");
	// 	fp=fopen(filename, "w");
	// 	fprintf(fp, "# T  specific_heat\n");
	// 	for(int it=1; it<N_T; it++){
	// 		double T_dif = temperature[it] - temperature[it-1];
	// 		double heat = (T_energy[it] - T_energy[it-1]) / T_dif;
	// 		double heat_err = sqrt( pow(T_energy_err[it], 2) + pow(T_energy_err[it-1], 2) ) / T_dif;
	// 		
	// 		fprintf(fp, "%.5e %.5e %.5e\n", temperature[it] - 0.5*T_dif, heat, heat_err);
	// 	}
	// 	fclose(fp);
	// 	printf("\n '%s'\n", filename);
		
		
		clock_t time_end = clock();
		print_time(time_start, time_end);
	// 	print_time(time_start, time_end, LOG_FILE);
	}
	
}



main(int argc, char* argv[])
{
	#if HYB_QMC_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	#endif // HYB_QMC_MPI
	
	hybqmc_init(&PQ, &SP, &TP_tau, &TP, &TP_sp, &TP_ch, &TP_tr, &D);
	
	init_params();
	
	hybqmc_set_nmc(n_mc);
	
	G = new dmft_green_func [N_S];
	
	switch(i_program){
		case 0:
			main_sub(0);
			break;
		case 1:
			main_T();
			break;
// 		case 2:
// 			main_H();
// 			break;
	}
	
	delete G;
	
	hybqmc_final();
	
	#if HYB_QMC_MPI
	MPI_Finalize();
	#endif // HYB_QMC_MPI
}
