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

int my_rank=0;

// ============================================================================

class InputParams{
public:
	// [model]
	int n_s;
	string file_Delta_tau, file_Delta_iw;  // mutually exclusive
	string file_Vsq;  // only when file_Delta_iw is set
	string file_U, file_ef;
	double beta;

	// [control]
	int n_tau;
	int n_tp, n_tp2;
	int rand_seed;
	int max_order;
	bool flag_tp;

	// [MC]
	int n_bin, n_msr, n_add, n_shift, n_warmup;
	double r_add, r_shift;

	// methods
	InputParams(bool verbose) : verbose(verbose) {};
	void read_params(string& file_ini);
	void summary();

private:
	bool verbose;
};

void InputParams::read_params(string& file_ini)
{
    boost::property_tree::ptree pt;

	if(verbose){
		cout << "\nRead file '" << file_ini << "'" << endl;
	}

	try{
	    boost::property_tree::read_ini(file_ini, pt);

		// [model]
		n_s = pt.get<int>("model.n_s");
		file_U = pt.get<string>("model.file_U");
		file_ef = pt.get<string>("model.file_ef");
		beta = pt.get<double>("model.beta");
		file_Delta_tau = pt.get<string>("model.file_Delta_tau", "");
		file_Delta_iw = pt.get<string>("model.file_Delta_iw", "");

		// file_Delta_tau and file_Delta_iw are mutually exclusive
		if (file_Delta_tau.empty() == file_Delta_iw.empty()){  // not XOR
			cerr << "INPUT_ERROR: 'file_Delta_tau' and 'file_Delta_iw' are mutuallly exclusive required." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (!file_Delta_iw.empty()){
			file_Vsq = pt.get<string>("model.file_Vsq");
		}

		// [control]
		n_tau = pt.get<int>("control.n_tau");
		rand_seed = pt.get<int>("control.rand_seed", 0);
		max_order = pt.get<int>("control.max_order", 1024);
		flag_tp = pt.get<bool>("control.flag_tp", false);
		n_tp = pt.get<int>("control.n_tp", 32);
		n_tp2 = pt.get<int>("control.n_tp2", 256);

		// [MC]
		n_warmup = pt.get<int>("MC.n_warmup", 1000000);
		n_bin = pt.get<int>("MC.n_bin", 10);
		n_msr = pt.get<int>("MC.n_msr");
		n_add = pt.get<int>("MC.n_add", 0);
		r_add = pt.get<double>("MC.r_add", 0.4);
		n_shift = pt.get<int>("MC.n_shift", 0);
		r_shift = pt.get<double>("MC.r_shift", 0.4);
	}
	catch(boost::property_tree::ptree_error& e){
		cerr << "INPUT_ERROR: " << e.what() << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

void InputParams::summary()
{
	if(verbose){
		cout << boolalpha;  // print bool type as true or false
		cout << "\nSummary of input parameters" << endl;
		cout << "=====================================" << endl;
		cout << "[model]" << endl;
		cout << "n_s = " << n_s << endl;
		cout << "file_Delta_tau = " << file_Delta_tau << endl;
		cout << "file_Delta_iw = " << file_Delta_iw << endl;
		cout << "file_Vsq = " << file_Vsq << endl;
		cout << "file_U = " << file_U << endl;
		cout << "file_ef = " << file_ef << endl;
		cout << "beta = " << beta << endl;

		cout << "\n[control]" << endl;
		cout << "n_tau = " << n_tau << endl;
		cout << "n_tp = " << n_tp << endl;
		cout << "n_tp2 = " << n_tp2 << endl;
		cout << "rand_seed = " << rand_seed << endl;
		cout << "max_order = " << max_order << endl;
		cout << "flag_tp = " << flag_tp << endl;

		cout << "\n[MC]" << endl;
		cout << "n_warmup = " << n_warmup << endl;
		cout << "n_bin = " << n_bin << endl;
		cout << "n_msr = " << n_msr << endl;
		cout << "n_add = " << n_add << endl;
		cout << "r_add = " << r_add << endl;
		cout << "n_shift = " << n_shift << endl;
		cout << "r_shift = " << r_shift << endl;
		cout << "=====================================" << endl;
	}
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

vec_vec_d load(const string& file_data, int n_col)
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

	// Count the number of rows
	int n_row=0;
	double data_line[256];
	while( read_data(fp, data_line) != 0 ){
		++n_row;
	}
	rewind(fp);

	vec_vec_d data;
	resize(data, n_row, n_col);

	for(int i=0; i<n_row; i++){
		if( read_data(fp, data_line) != n_col ){
			// printf("ERROR: Invalid file '%s'\n", file_data);
			cout << "ERROR: Invalid data structure '" << file_data << "'" << endl;
			fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for(int j=0; j<n_col; j++){
			data[i][j] = data_line[j];
		}
	}
	return data;
}

void print_size(const vec_vec_d& data)
{
	cout << "  (" << data.size() << ", " << data[0].size() << ") loaded." << endl;
}

void verify_size(const vec_vec_d& data, size_t n_rows, size_t n_cols)
{
	if( !check_size(data, n_rows, n_cols) ){
		cerr << "Invalid data size : Require (" << n_rows << ", " << n_cols << ")" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

vec_d read_ef(const string& file_data, int n_s)
{
	vec_vec_d data = load(file_data, 1);

	if(my_rank==0){
		cout << "\nRead file '" << file_data << "'" << endl;
		print_size(data);
		verify_size(data, n_s, 1);
	}

	vec_d ef(n_s);
	for(int i=0; i<n_s; i++){
		ef[i] = data[i][0];
	}
	return ef;
}

vec_vec_d read_U(const string& file_data, int n_s)
{
	vec_vec_d data = load(file_data, n_s);

	if(my_rank==0){
		cout << "\nRead file '" << file_data << "'" << endl;
		print_size(data);
		verify_size(data, n_s, n_s);
	}

	return data;
}

vec_vec_d read_Delta_tau(const string& file_data, int n_s, vec_d& tau)
{
	vec_vec_d data = load(file_data, n_s+1);

	if(my_rank==0){
		cout << "\nRead file '" << file_data << "'" << endl;
		print_size(data);
	}

	int n_t = data.size() - 1;

	tau.resize(n_t + 1);
	for(int i=0; i<=n_t; i++){
		tau[i] = data[i][0];
	}

	vec_vec_d Delta_tau;
	resize(Delta_tau, n_s, n_t + 1);
	for(int s=0; s<n_s; s++){
		for(int i=0; i<=n_t; i++){
			Delta_tau[s][i] = data[i][1+s];
		}
	}
	return Delta_tau;
}

vec_vec_c read_Delta_iw(const string& file_data, int n_s)
{
	vec_vec_d data = load(file_data, 2*n_s);

	if(my_rank==0){
		cout << "\nRead file '" << file_data << "'" << endl;
		print_size(data);
	}

	int n_w = data.size();

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
	vec_vec_d data = load(file_data, 1);

	if(my_rank==0){
		cout << "\nRead file '" << file_data << "'" << endl;
		print_size(data);
		verify_size(data, n_s, 1);
	}

	vec_d Vsq(n_s);
	for(int i=0; i<n_s; i++){
		Vsq[i] = data[i][0];
	}
	return Vsq;
}


// ============================================================================

void print_pq(hyb_qmc_params& prm, phys_quant& PQ, t_sp& SP)
{
	int N_S = SP.size();

	FILE *fp;
	char filename[128];

	sprintf(filename, "x.dat");
	fp=fopen(filename, "w");

	fprintf(fp, "%.5e", 1./prm.beta);
	for(int s=0; s<N_S; s++)  fprintf(fp, " %.5e %.3e", SP[s].f_number, SP[s].f_number_err);
	fprintf(fp, " %.5e %.4e", PQ.stat_suscep_sp, PQ.stat_suscep_sp_err);
	fprintf(fp, " %.5e %.4e", PQ.stat_suscep_ch, PQ.stat_suscep_ch_err);
	// fprintf(fp, " %.5e %.4e", TP[0][1].chi_tau[0], TP[0][1].chi_tau_err[0]);
	for(int s=0; s<N_S; s++){
		fprintf(fp, " %.5e %.5e", real(SP[s].self_omega_dyson[0]), imag(SP[s].self_omega_dyson[0]));
		fprintf(fp, " %.5e %.5e", real(SP[s].self_omega_dyson[1]), imag(SP[s].self_omega_dyson[1]));
	}
	fprintf(fp, "\n");

	fclose(fp);
	printf("\n '%s'\n", filename);
}

void print_statistics(phys_quant& PQ)
{
	FILE *fp;
	char filename[128];

	// Find the maximum expansion order
	int max_k = 0;
	for(int i=0; i<PQ.Z_k[0].size(); i++){
		for(int s=0; s<PQ.Z_k.size(); s++){
			if( PQ.Z_k[s][i] != 0. )  max_k = i;
		}
	}

	sprintf(filename, "stat.dat");
	fp=fopen(filename, "w");
	// for(int i=0; i<PQ.Z_k[0].size(); i++){
	for(int i=0; i<=max_k; i++){
		fprintf(fp, "%d", i);
		for(int s=0; s<PQ.Z_k.size(); s++){
			if( PQ.Z_k[s][i] != 0. )  fprintf(fp, " %.5e", PQ.Z_k[s][i]);
			else  fprintf(fp, " ?");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
}

void print_single_particle(hyb_qmc_params& prm, phys_quant& PQ, t_sp& SP)
{
	int N_S = SP.size();
	int N_TAU = SP[0].Gf_tau.size() - 1;

	printf("\n Occupation number           (measured at tau=0)\n");
	for(int s=0; s<N_S; s++){
		printf("    %d: %lf +- %lf", s, SP[s].f_number, SP[s].f_number_err);
		printf("    (%lf +- %lf)\n", SP[s].f_number2, SP[s].f_number2_err);
	}
	printf("  tot: %lf +- %lf\n", PQ.occup_tot, PQ.occup_tot_err);

	FILE *fp;
	char filename[128];

	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}

	sprintf(filename, "Gf_t.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.4e", tau[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", SP[s].Gf_tau[i], SP[s].Gf_tau_err[i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);

	double iw[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		iw[i] = (double)(2*i+1) * M_PI / prm.beta;
	}

	sprintf(filename, "Gf_w.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%.4e", iw[i]);
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
	// 		pade_init(i_omega_f, SP[s].Gf_w.data(), Gf_pade[s], N_PADE);
	// 	}

	// 	sprintf(filename, "Gf_pade.dat");
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


	sprintf(filename, "self_w.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
// 	for(int i=0; i<100; i++){
		fprintf(fp, "%.4e", iw[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_omega[i]), imag(SP[s].self_omega[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);


	sprintf(filename, "self_w_dyson.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
// 	for(int i=0; i<100; i++){
		fprintf(fp, "%.4e", iw[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(SP[s].self_omega_dyson[i]), imag(SP[s].self_omega_dyson[i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);


	sprintf(filename, "GSigma_t.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.4e", tau[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", SP[s].GSigma_tau[i], SP[s].GSigma_tau_err[i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);
}


void print_two_particle(phys_quant& PQ, t_tp& TP, vec_d& TP_tau, two_particle& TP_sp, two_particle& TP_ch)
{
	printf("\n static susceptibility\n");
	printf("  spin   : %.5e +- %.4e\n", PQ.stat_suscep_sp, PQ.stat_suscep_sp_err);
	printf("  charge : %.5e +- %.4e\n", PQ.stat_suscep_ch, PQ.stat_suscep_ch_err);

	int N_S = TP.size();
	int N_TP = TP[0][0].chi_tau.size() - 1;
	int N_TP2 = TP[0][0].chi_omega.size() - 1;

	FILE *fp;
	char filename[128];

	sprintf(filename, "chi_t.dat");
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
	printf("\n '%s'\n", filename);



	sprintf(filename, "chi_w.dat");
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

	// sprintf(filename, "chi_pade.dat");
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


void print_Delta_tau(hyb_qmc_params& prm, vec_vec_d& delta_tau, vec_d& tau)
{
	int N_S = delta_tau.size();
	int N_TAU = delta_tau[0].size() - 1;

	FILE *fp;
	char filename[128];

	sprintf(filename, "delta_t.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.4e", tau[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e", delta_tau[s][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("\n '%s'\n", filename);

	printf(" Time resolution  : %.2e\n", tau[1] - tau[0]);
}


void print_Delta_iw(hyb_qmc_params& prm, vec_vec_c& delta_omega, vec_d& Vsq)
{
	int N_S = delta_omega.size();
	int N_TAU = delta_omega[0].size() * 2;

	FILE *fp;
	char filename[128];

	double iw[N_TAU/2];
	for(int i=0; i<N_TAU/2; i++){
		iw[i] = (double)(2*i+1) * M_PI / prm.beta;
	}

	double tau[N_TAU+1];
	for(int i=0; i<=N_TAU; i++){
		tau[i] = (double)i * prm.beta / (double)N_TAU;
	}

	printf("\n Frequency cutoff : %.2e\n", iw[N_TAU/2-1]);
	printf(" Time resolution  : %.2e\n", tau[1] - tau[0]);

	sprintf(filename, "delta_w.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<N_TAU/2; i++){
		fprintf(fp, "%.4e", iw[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e %.6e", real(delta_omega[s][i]), imag(delta_omega[s][i]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf(" '%s'\n", filename);


	// FFT: G(iw) -> G(tau)
	vec_vec_d delta_tau;
	delta_tau.resize(N_S);
	for(int s=0; s<N_S; s++){
		fft_fermion_iw2tau(delta_tau[s], delta_omega[s], prm.beta, Vsq[s]);
		assert (delta_tau[s].size() == N_TAU);
		delta_tau[s].push_back(- Vsq[s] - delta_tau[s][0]);  // [N_TAU]
	}

	sprintf(filename, "delta_t.dat");
	fp=fopen(filename, "w");
	for(int i=0; i<=N_TAU; i++){
		fprintf(fp, "%.4e", tau[i]);
		for(int s=0; s<N_S; s++){
			fprintf(fp, " %.6e", delta_tau[s][i]);
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

	if(my_rank==0){
		cout << "HYBQMC start" << endl;
	}
	clock_t time_start = clock();

	// Input parameters
	InputParams in(my_rank==0);
	in.read_params(file_ini);
	in.summary();

	// Initialize QMC worker
	HybQMC MC(in.max_order, in.n_s, in.n_tau, in.n_tp, in.n_tp2, in.rand_seed);

	// Set MC parameters
	num_mc n_mc;
	n_mc.N_MSR = in.n_msr;
	n_mc.N_BIN = in.n_bin;
	n_mc.N_ADD = in.n_add;
	n_mc.N_SHIFT = in.n_shift;
	n_mc.N_WARMUP = in.n_warmup;
	n_mc.R_ADD = in.r_add;
	n_mc.R_SHIFT = in.r_shift;
	MC.set_nmc(n_mc);

	// Set parameters
	hyb_qmc_params prm;
	prm.beta = in.beta;
	prm.ef = read_ef(in.file_ef, in.n_s);
	prm.U = read_U(in.file_U, in.n_s);
	MC.set_params(prm);

	// Read and set Delta(iw)
	if(!in.file_Delta_iw.empty()){
		vec_vec_c Delta_omega = read_Delta_iw(in.file_Delta_iw, in.n_s);
		vec_d Vsq = read_Vsq(in.file_Vsq, in.n_s);
		if(my_rank==0){
			print_Delta_iw(prm, Delta_omega, Vsq);
		}
		MC.set_Delta_iw(Delta_omega, Vsq);
	}
	// Read and set Delta(tau)
	if(!in.file_Delta_tau.empty()){
		vec_d tau;
		vec_vec_d Delta_tau = read_Delta_tau(in.file_Delta_tau, in.n_s, tau);
		if(my_rank==0){
			print_Delta_tau(prm, Delta_tau, tau);
		}
		MC.set_Delta_tau(Delta_tau, tau);
	}

	// QMC calc
	int flag_tp = in.flag_tp;
	MC.eval(flag_tp);

	// Get results
	t_sp SP = MC.get_SP();
	t_tp TP = MC.get_TP();
	vec_d TP_tau = MC.get_TP_tau();
	two_particle TP_sp = MC.get_TP_sp();
	two_particle TP_ch = MC.get_TP_ch();
	phys_quant PQ = MC.get_PQ();

	// Output results
	if(my_rank==0){
		print_pq(prm, PQ, SP);
		print_statistics(PQ);
		print_single_particle(prm, PQ, SP);

		if(flag_tp){
			print_two_particle(PQ, TP, TP_tau, TP_sp, TP_ch);
		}

		clock_t time_end = clock();
		cout << "\nHYBQMC finish" << endl;
		char str[100];
		sprint_time(str, time_end - time_start);
		printf("  time:%s", str);
	}

	#if HYB_QMC_MPI
	MPI_Finalize();
	#endif // HYB_QMC_MPI
}
