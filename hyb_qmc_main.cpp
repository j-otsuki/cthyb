/*

Continuous-time quantum Monte Carlo
for the impurity Anderson model (hybridization expansion)

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

static const char tag[64] = "v1.12";

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/optional.hpp>
#include "hyb_qmc.h"
#include "vector_type.hpp"
// #include "green_func_0.h"
// #include "dmft.h"
#include "fft.h"
#include "pade.h"

using namespace std;

// #define DMFT_ISO 1
// 0: inisotropic bath, 1: isotropic bath

// #define DMFT_OUTPUT_ITER 0
// output cavity Green function at each iteraction

// int DMFT_ITER=0;
// double DMFT_UPDATE_RATE=1.0;

// int N_PADE = std::min(N_TAU/2, 8192);


// struct hyb_qmc_params prm;
// double prm_ave_n, prm_chem_pot;
// double prm_V_sqr[N_S], prm_V_sqr_imp[N_S], prm_D, prm_f_disp=0;
// double prm_E_f[N_S], prm_E_c[N_S];

// struct num_mc n_mc;

// struct phys_quant *PQ;
// struct single_particle *SP;
// struct two_particle **TP, *TP_sp, *TP_ch;
// struct two_particle_tr *TP_tr;
// double *TP_tau;
// struct phonon_g *D;

// complex<double> G0_omega[N_S][N_TAU/2];

// struct dmft_green_func G;
// struct dmft_green_func *G;

// double SP_c_number[N_S];

int my_rank;


int i_program;
// char str_program[3][20]={
// 	"T,H fixed",
// 	"T varies",
// 	"H varies"
// };

int flag_tp;
char str_flag[2][100]={
	"single-particle Green function",  // 0
	"single-particle & two-particle Green function",  // 1
};

// int flag_dmft;
// char str_dmft[2][100]={
// 	"impurity problem",
// 	"DMFT",
// };

// int flag_ensemble;
// char str_ensemble[2][100]={
// 	"grand canonical ensemble (chemical potential given)",
// 	"canonical ensemble (band filling given)"
// };



// ============================================================================

/*
void init_params()
{
	FILE *fp;
	if( (fp = fopen("hyb_qmc_param.init", "r")) == NULL ){
		printf("'hyb_qmc_param.init' not found\n");
		exit(0);
	}

	if(my_rank==0){
		printf("\nCT-QMC for the Anderson model (CT-HYB)  %s\n", tag);
		// printf(" DOS = %d\n", DOS);
	}

	double data[16];

	if( read_data(fp, data) != 1 )  exit(0);
	i_program = (int)data[0];

	if( read_data(fp, data) != 1 )  exit(0);
	flag_tp = (int)data[0];

	switch( read_data(fp, data) ){
		// case 4:
		// 	prm_f_disp = data[3];
		// case 3:
		// 	DMFT_UPDATE_RATE = data[2];
		// case 2:
		// 	DMFT_ITER = (int)data[1];
		// case 1:
		// 	flag_dmft = (int)data[0];
		// 	break;
		default:
			// exit(0);
			break;
	}

	if( read_data(fp, data) != 2 )  exit(0);
	// flag_ensemble = (int)data[0];
	// if( flag_ensemble==0 ){  // mu given
	// 	prm_chem_pot = data[1];
	// }
	// else{  // n given
	// 	prm_ave_n = data[1];
	// 	prm_chem_pot = 0;
	// }

	double prm_epsilon_f, prm_H;

	double prm_U;
	if( fscanf(fp, "%lf %lf %lf %lf", &prm_D, &prm_V_sqr[0], &prm_epsilon_f, &prm_U) != 4 )  exit(0);
	if( fscanf(fp, "%lf %lf", &prm.beta, &prm_H) != 2 )  exit(0);
	if( isinf(prm_U) )  prm.UINF = 1;

	// if(PHONON){
	// 	if( fscanf(fp, "%lf %lf", &prm.w_ph, &prm.g) != 2 )  exit(0);
	// }

	if(my_rank==0){
		// printf("\n %d: %s\n", i_program, str_program[i_program]);
		printf("\n %d:\n", i_program);
		printf(" %d: %s\n", flag_tp, str_flag[flag_tp]);
		// printf(" %d: %s\n", flag_dmft, str_dmft[flag_dmft]);
		// printf(" %d: %s\n", flag_ensemble, str_ensemble[flag_ensemble]);

		printf("\n D = %.4e\n", prm_D);
		printf(" V_sqr = %.4e\n", prm_V_sqr[0]);
		printf(" epsilon_f = %.4e\n", prm_epsilon_f);
		printf(" U = %.4e  (UINF = %d)\n", prm_U, prm.UINF);
		printf(" beta = %.4e\n", prm.beta);
		printf(" H = %.4e\n", prm_H);
		// if(PHONON){
		// 	printf(" w_ph = %.4e\n", prm.w_ph);
		// 	printf(" g = %.4e\n", prm.g);
		// }

		// if( !flag_ensemble )  printf(" chem_pot = %.4e\n", prm_chem_pot);  // grand canonical
		// else                  printf(" ave_nc = %.4e\n", prm_ave_n);  // canonical

		// if( flag_dmft==1 ){
		// 	printf(" DMFT_ITER = %d\n", DMFT_ITER);
		// 	printf(" DMFT_UPDATE_RATE = %.2lf\n", DMFT_UPDATE_RATE);
		// 	printf(" f_disp = %.3lf\n", prm_f_disp);
		// }
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

	// if(i_program==1){
	// 	if( fscanf(fp, "%d %d %lf %lf", &i_tmesh, &N_T, &T_MIN, &T_MAX) != 4 )  exit(0);
	// }

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
*/


class InputParams{
public:
	// [model]
	int n_s;
	string file_Delta, file_U, file_ef, file_Vsq;
	double beta;

	// [control]
	int n_tau;
	int n_tp, n_tp2;
	int rand_seed;
	int max_order;

	// [MC]
	int n_bin, n_msr, n_add, n_shift;

	// methods
	void read_params(string& file_ini);
	void summary();
};

void InputParams::read_params(string& file_ini)
{
    boost::property_tree::ptree pt;

	cout << "Reading file '" << file_ini << "'" << endl;
	try{
	    boost::property_tree::read_ini(file_ini, pt);

		// [model]
		n_s = pt.get<int>("model.n_s");
		file_Delta = pt.get<string>("model.file_Delta");
		file_Vsq = pt.get<string>("model.file_Vsq");
		file_U = pt.get<string>("model.file_U");
		file_ef = pt.get<string>("model.file_ef");
		beta = pt.get<double>("model.beta");

		// [control]
		n_tau = pt.get<int>("control.n_tau");
		n_tp = pt.get<int>("control.n_tp");
		n_tp2 = pt.get<int>("control.n_tp2");
		rand_seed = pt.get<int>("control.rand_seed", 0);
		max_order = pt.get<int>("control.max_order", 1024);

		// [MC]
		n_bin = pt.get<int>("MC.n_bin", 10);
		n_msr = pt.get<int>("MC.n_msr");
		n_add = pt.get<int>("MC.n_add", -4);
		n_shift = pt.get<int>("MC.n_shift", -4);
	}
	catch(boost::property_tree::ptree_error& e){
		cerr << "ERROR: " << e.what() << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

void InputParams::summary()
{
	cout << "\nSummary of input parameters" << endl;
	cout << "=====================================" << endl;
	cout << "[model]" << endl;
	cout << "n_s = " << n_s << endl;
	cout << "file_Delta = " << file_Delta << endl;
	cout << "file_Vsq = " << file_Vsq << endl;
	cout << "file_U = " << file_U << endl;
	cout << "file_ef = " << file_ef << endl;
	cout << "beta = " << beta << endl;

	cout << "\n[control]" << endl;
	cout << "n_tau = " << n_tau << endl;
	cout << "n_tp = " << n_tp << endl;
	cout << "n_tp2 = " << n_tp2 << endl;
	cout << "rand_seed = " << rand_seed << endl;
	cout << "max_order = " << rand_seed << endl;

	cout << "\n[MC]" << endl;
	cout << "n_bin = " << n_bin << endl;
	cout << "n_msr = " << n_msr << endl;
	cout << "n_add = " << n_add << endl;
	cout << "n_shift = " << n_shift << endl;
	cout << "=====================================" << endl;
}

// ============================================================================

// return data[row][col]
// vec_vec_d read_data(string& file_data, int n_row, int n_col)
// {
// 	vec_vec_d data(n_row, n_col);
// 	ifstream ifs(file_data);
// 	if (!ifs) {
// 		cerr << "ERROR: File '"<< file_data << "'not found" << endl;
// 		MPI_Abort(MPI_COMM_WORLD, 1);
// 	}

// 	for(int i=0; i<n_row; i++){
// 		std::string buf;
// 		getline(ifs, buf);
// 		cout << buf << endl;
// 		for(int j=0; j<n_col; j++){
// 			sscanf("%lf", &data[i][j]);
// 		}
// 	}
// }
vec_vec_d load(const string& file_data, int n_row, int n_col)
{
	FILE *fp;
	// if( (fp = fopen(file_data.c_str(), "r")) == NULL ){
	// 	printf("'%s' not found\n", file_data.c_str());
	// 	MPI_Abort(MPI_COMM_WORLD, 1);
	// }
	if( (fp = fopen(file_data.c_str(), "r")) == NULL ){
		// printf("'%s' not found\n", file_data);
		cout << "'" << file_data << "' not found" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	vec_vec_d data;
	resize(data, n_row, n_col);

	double data_line[256];
	for(int i=0; i<n_row; i++){
		if( read_data(fp, data_line) != n_col ){
			// printf("ERROR: Invalid file '%s'\n", file_data);
			cout << "ERROR: Invalid file '" << file_data << "'" << endl;
			fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for(int j=0; j<n_col; j++){
			data[i][j] = data_line[j];
		}
	}
	return data;
}

vec_d read_ef(const string& file_data, int n_s)
{
	cout << "\nRead file '" << file_data << "'" << endl;
	vec_vec_d data = load(file_data, n_s, 1);
	cout << "  done" << endl;

	vec_d ef(n_s);
	for(int i=0; i<n_s; i++){
		ef[i] = data[i][0];
	}
	return ef;
}

vec_vec_d read_U(const string& file_data, int n_s)
{
	cout << "\nRead file '" << file_data << "'" << endl;
	vec_vec_d data = load(file_data, n_s, n_s);
	cout << "  done" << endl;

	return data;
}

vec_vec_c read_Delta(const string& file_data, int n_s, int n_w)
{
	cout << "\nRead file '" << file_data << "'" << endl;
	// cout << "n_s=" << n_s << " , n_w=" << n_w << endl;
	vec_vec_d data = load(file_data, n_w, 2*n_s);
	cout << "  done" << endl;

	vec_vec_c Delta_omega;
	resize(Delta_omega, n_s, n_w);
	for(int s=0; s<n_s; s++){
		for(int i=0; i<n_w; i++){
			Delta_omega[s][i] = complex<double>(data[i][2*s], data[i][2*s+1]);
		}
	}
	return Delta_omega;
}

vec_d read_Vsq(const string& file_data, int n_s)
{
	cout << "\nRead file '" << file_data << "'" << endl;
	vec_vec_d data = load(file_data, n_s, 1);
	cout << "  done" << endl;

	vec_d Vsq(n_s);
	for(int i=0; i<n_s; i++){
		Vsq[i] = data[i][0];
	}
	return Vsq;
}


// ============================================================================

void print_pq(int num, hyb_qmc_params& prm, phys_quant& PQ, t_sp& SP)
{
	int N_S = SP.size();

	char filename[128];
	sprintf(filename, DATA_DIR "T_depend.dat");

	static FILE *fp=fopen(filename, "w");
	fp=fopen(filename, "a");

	fprintf(fp, "%.5e", 1./prm.beta);
	for(int s=0; s<N_S; s++)  fprintf(fp, " %.5e %.3e", SP[s].f_number, SP[s].f_number_err);
	fprintf(fp, " %.5e %.4e", PQ.stat_suscep_sp, PQ.stat_suscep_sp_err);
	fprintf(fp, " %.5e %.4e", PQ.stat_suscep_ch, PQ.stat_suscep_ch_err);
	// fprintf(fp, " %.5e %.4e", TP[0][1].chi_tau[0], TP[0][1].chi_tau_err[0]);
	for(int s=0; s<N_S; s++){
		fprintf(fp, " %.5e %.5e", real(SP[s].self_f[0]), imag(SP[s].self_f[0]));
		fprintf(fp, " %.5e %.5e", real(SP[s].self_f[1]), imag(SP[s].self_f[1]));
	}
	// for(int s=0; s<N_S; s++)  fprintf(fp, " %.5e", SP_c_number[s]);
// 	if( flag_dmft ){
// 		fprintf(fp, " %.6e %.6e", prm_chem_pot, prm_ave_n);
// 		for(int s=0; s<N_S; s++){
// 			fprintf(fp, " %.5e %.5e", G[s].n_f, G[s].n_c);
// // 			fprintf(fp, " %.5e %.5e", real(G[s].F[0])*prm_V_sqr[s], imag(G[s].F[0])*prm_V_sqr[s]);
// // 			fprintf(fp, " %.5e %.5e", real(G[s].F[1])*prm_V_sqr[s], imag(G[s].F[1])*prm_V_sqr[s]);
// 			fprintf(fp, " %.5e %.5e", real(G[s].self_c[0]), imag(G[s].self_c[0]));
// 			fprintf(fp, " %.5e %.5e", real(G[s].self_c[1]), imag(G[s].self_c[1]));
// 		}
// 	}
// 	if( PHONON ){
// 		fprintf(fp, " %.6e", D->occup);
// 	}

	fprintf(fp, "\n");
	fclose(fp);
	printf("\n '%s' updated\n", filename);
}

void print_statistics(int num, phys_quant& PQ)
{
	FILE *fp;
	char filename[40];

	sprintf(filename, DATA_DIR "%02d-stat.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<PQ.Z_k[0].size(); i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<PQ.Z_k.size(); s++){
			if( PQ.Z_k[s][i] != 0. )  fprintf(fp, " %.5e", PQ.Z_k[s][i]);
			else  fprintf(fp, " ?");
		}
		fprintf(fp, "\n");
	}
	printf("\n '%s'\n", filename);
	fclose(fp);

	sprintf(filename, DATA_DIR "%02d-stat_tot.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<PQ.Z_ktot.size(); i++){
		fprintf(fp, "%d", i);
		if( PQ.Z_ktot[i] != 0. )  fprintf(fp, " %.5e", PQ.Z_ktot[i]);
		else  fprintf(fp, " ?");
		fprintf(fp, "\n");
	}
	printf("\n '%s'\n", filename);
	fclose(fp);
}

void print_single_particle(int num, hyb_qmc_params& prm, phys_quant& PQ, t_sp& SP)
{
	int N_S = SP.size();
	int N_TAU = SP[0].Gf_tau.size() - 1;

	printf("\n f number                   (measured at tau=0)\n");
	for(int s=0; s<N_S; s++){
		printf("   %d: %lf +- %lf", s, SP[s].f_number, SP[s].f_number_err);
		printf("   (%lf +- %lf)\n", SP[s].f_number2, SP[s].f_number2_err);
	}
	printf(" tot: %lf +- %lf\n", PQ.occup_tot, PQ.occup_tot_err);

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


	// Number of frequencies used for Pade approximation of Gf (N_PADE^2 memory)
	// MAX: N_TAU/2 (all data),  MIN: 0 (do not evaluate)
	// int N_PADE = std::min(N_TAU/2, 8192);

	// complex<double> i_omega_f[N_PADE];

	// if(N_PADE){
	// 	complex<double> Gf_pade[N_S][N_PADE];

	// 	for(int i=0; i<N_PADE; i++){
	// 		double omega_f = (double)(2*i+1) * M_PI / prm.beta;
	// 		i_omega_f[i] = IMAG * omega_f;
	// 	}

	// 	for(int s=0; s<N_S; s++){
	// 		pade_init(i_omega_f, SP[s].Gf_omega.data(), Gf_pade[s], N_PADE);
	// 	}

	// 	sprintf(filename, DATA_DIR "%02d-Gf_pade.dat", num);
	// 	fp=fopen(filename, "w");
	// 	for(int i=0; i<1000; i++){

	// 		complex<double> w = -2 * prm_D + 4 * prm_D / 999. * (double)i;
	// 		copmlex<double> w =
	// 		// double delta0 = M_PI * prm_V_sqr_imp[0] / prm_D * 0.5;

	// 		fprintf(fp, "%.6e", real(w));
	// 		for(int s=0; s<N_S; s++){
	// 			complex<double> temp_Gf = pade_calc(i_omega_f, Gf_pade[s], w, N_PADE);
	// 			fprintf(fp, " %.6e %.6e", real(temp_Gf), imag(temp_Gf));
	// 		}
	// 		// fprintf(fp, " %.6e",
	// 		//  -delta0 / (pow(real(w)-prm.ef[0]-prm.U[0][0]*0.5,2) + pow(delta0,2)) );

	// 		fprintf(fp, "\n");
	// 	}
	// 	fclose(fp);
	// 	printf(" '%s'\n", filename);
	// }


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

}


void print_two_particle(int num, phys_quant& PQ, t_tp& TP, vec_d& TP_tau, two_particle& TP_sp, two_particle& TP_ch)
{
	printf("\n static susceptibility\n");
	printf("  spin   : %.5e +- %.4e\n", PQ.stat_suscep_sp, PQ.stat_suscep_sp_err);
	printf("  charge : %.5e +- %.4e\n", PQ.stat_suscep_ch, PQ.stat_suscep_ch_err);

// 	printf("\n static susceptibility (test)\n");
// 	printf("  spin   : %.5e\n", real(TP_sp.chi_omega[0]));
// 	printf("  charge : %.5e\n", real(TP_ch.chi_omega[0]));


	int N_S = TP.size();
	int N_TP = TP[0][0].chi_tau.size() - 1;
	int N_TP2 = TP[0][0].chi_omega.size() - 1;

	FILE *fp;
	char filename[40];

// 	double tau[N_TAU+1];
// 	for(int i=0; i<=N_TAU; i++){
// 		tau[i] = (double)i * prm.beta / (double)N_TAU;
// 	}

	//
	// non-interacting
	//
	// double Gf0_tau[N_TAU+1];
	// Gf_calc_fft(Gf0_tau, prm.beta, prm_D, prm.ef[0], prm.U[0][0], prm_V_sqr_imp[0]);
	// Gf0_tau[N_TAU] = -1.0 - Gf0_tau[0];

	// double chi0_tau[N_TAU+1];
	// for(int i=0; i<=N_TAU/2; i++)  chi0_tau[i] = Gf0_tau[i] * Gf0_tau[N_TAU-i];

	// sprintf(filename, DATA_DIR "%02d-chi0.dat", num);
	// fp=fopen(filename, "w");
	// for(int i=0; i<=N_TAU/2; i++){
	// 	double tau = (double)i * prm.beta / (double)N_TAU;
	// 	fprintf(fp, "%.5e %.5e\n", tau, chi0_tau[i]);
	// }
	// fclose(fp);
	// printf("\n '%s'\n", filename);


	sprintf(filename, DATA_DIR "%02d-chi_tau.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TP; i++){
// 		double tau = (double)i * prm.beta / (double)(2*N_TP);
		fprintf(fp, "%.5e", TP_tau[i]);
		fprintf(fp, " %.5e %.4e", TP_sp.chi_tau[i], TP_sp.chi_tau_err[i]);
		fprintf(fp, " %.5e %.4e", TP_ch.chi_tau[i], TP_ch.chi_tau_err[i]);
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
		fprintf(fp, "%d %.5e %.5e\n", i, real(TP_sp.chi_omega[i]), real(TP_ch.chi_omega[i]));
	}
	fclose(fp);
	printf(" '%s'\n", filename);


	// complex<double> chi_sp_pade[N_TP2], chi_ch_pade[N_TP2], i_omega_b[N_TP2];

	// for(int i=0; i<N_TP2; i++)  i_omega_b[i] = IMAG * (double)(2*i) * M_PI / prm.beta;

	// pade_init(i_omega_b, TP_sp.chi_omega.data(), chi_sp_pade, N_TP2);
	// pade_init(i_omega_b, TP_ch.chi_omega.data(), chi_ch_pade, N_TP2);

	// sprintf(filename, DATA_DIR "%02d-chi_pade.dat", num);
	// fp=fopen(filename, "w");
	// for(int i=0; i<1000; i++){
	// 	complex<double> w = 1.0 * prm_D / 999. * (double)i;

	// 	complex<double> temp_chi_sp = pade_calc(i_omega_b, chi_sp_pade, w, N_TP2);
	// 	complex<double> temp_chi_ch = pade_calc(i_omega_b, chi_ch_pade, w, N_TP2);

	// 	fprintf(fp, "%.6e %.6e %.6e %.6e %.6e\n", real(w),
	// 	 real(temp_chi_sp), imag(temp_chi_sp), real(temp_chi_ch), imag(temp_chi_ch));
	// }
	// fclose(fp);
	// printf(" '%s'\n", filename);


}


void print_Delta(int num, vec_vec_c& delta_omega)
{
	int N_S = delta_omega.size();
	int N_TAU = delta_omega[0].size() * 2;

	FILE *fp;
	char filename[40];

	// double tau[N_TAU+1];
	// for(int i=0; i<=N_TAU; i++){
	// 	tau[i] = (double)i * prm.beta / (double)N_TAU;
	// }

	// double G0_tau[N_TAU+1];
	// vec_vec_d delta_tau;
	// resize(delta_tau, N_S, N_TAU+1);
	// for(int s=0; s<N_S; s++){
	// 	fft_fermion_radix2_omega2tau(delta_tau[s], delta_omega[s], prm.beta, N_TAU, Vsq[s]);
	// }

	// sprintf(filename, DATA_DIR "%02d-delta_tau.dat", num);
	// fp=fopen(filename, "w");
	// for(int i=0; i<N_TAU; i++){
	// 	fprintf(fp, "%.4e %.6e\n", tau[i], G0_tau[i]);
	// 	for(int s=0; s<N_S; s++){
	// 		fprintf(fp, " %.6e", delta_tau[s][i]);
	// 	}
	// }
	// fclose(fp);
	// printf("\n '%s'\n", filename);


	sprintf(filename, DATA_DIR "%02d-delta_omega.dat", num);
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(delta_omega[s][i]), imag(delta_omega[s][i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);

}



// ============================================================================



int main(int argc, char* argv[])
{
	#if HYB_QMC_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	#endif // HYB_QMC_MPI

	string file_ini;
	switch(argc){
		case 2:
			file_ini = string(argv[1]);
			break;
		case 1:
			cerr << "ERROR: Missing argument. Input file name required." << endl;
			exit(1);
		default:
			cerr << "ERROR: Too many arguments. Input file name required." << endl;
			exit(1);
	}

	// init_params();
	InputParams in;
	// cout << file_ini << endl;
	// read_params(file_ini);
	in.read_params(file_ini);
	in.summary();

	// hybqmc_init(&PQ, &SP, &TP_tau, &TP, &TP_sp, &TP_ch, &TP_tr, &D);

	HybQMC Q(in.max_order, in.n_s, in.n_tau, in.n_tp, in.n_tp2, in.rand_seed);

	// hybqmc_set_nmc(n_mc);
	num_mc n_mc;
	n_mc.N_MSR = in.n_msr;
	n_mc.N_BIN = in.n_bin;
	n_mc.N_ADD = in.n_add;
	n_mc.N_SHIFT = in.n_shift;
	Q.set_nmc(n_mc);

	// G = new dmft_green_func [N_S];

	clock_t time_start = clock();

	hyb_qmc_params prm;
	prm.beta = in.beta;
	prm.ef = read_ef(in.file_ef, in.n_s);
	prm.U = read_U(in.file_U, in.n_s);
	Q.set_params(prm);

	// Delta_omega.resize(n_s, int(n_tau/2));
	// vec_vec_c Delta_omega;  // [N_S][N_TAU/2]
	// resize(Delta_omega, in.n_s, int(in.n_tau/2));
	// for(int s=0; s<N_S; s++)  G0_omega_calc(G0_omega[s], prm.beta, prm_D);
	//	G0_omega_calc(G0_omega, N_S, prm.beta, prm_D, flag_ensemble, prm_ave_n, prm_chem_pot, prm_E_c);
	// vector<double> prm_V_sqr(in.n_s);


	// hybqmc_set_G0(G0_omega, prm_V_sqr_imp);
	// read_Delta(in.file_Delta, in.n_s, int(in.n_tau/2));
	vec_vec_c Delta_omega = read_Delta(in.file_Delta, in.n_s, in.n_tau/2);
	vec_d Vsq = read_Vsq(in.file_Vsq, in.n_s);
	Q.set_Delta(Delta_omega, Vsq);

	int num = 0;
	if(my_rank==0){
		print_Delta(num, Delta_omega);
	}

	// hybqmc_eval(flag_tp);
	Q.eval(flag_tp);

	t_sp SP = Q.get_SP();
	t_tp TP = Q.get_TP();
	vec_d TP_tau = Q.get_TP_tau();
	two_particle TP_sp = Q.get_TP_sp();
	two_particle TP_ch = Q.get_TP_ch();
	phys_quant PQ = Q.get_PQ();

	if(my_rank==0){
		print_pq(num, prm, PQ, SP);
		print_statistics(num, PQ);
		print_single_particle(num, prm, PQ, SP);

		if(flag_tp){
			print_two_particle(num, PQ, TP, TP_tau, TP_sp, TP_ch);
		}

		clock_t time_end = clock();
		print_time(time_start, time_end);
	// 	print_time(time_start, time_end, LOG_FILE);
	}


	// delete G;

	// hybqmc_final();

	#if HYB_QMC_MPI
	MPI_Finalize();
	#endif // HYB_QMC_MPI
}
