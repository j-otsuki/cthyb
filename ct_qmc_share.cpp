/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "ct_qmc_share.h"


//=========================================================================================================

// // int [0:n)
// int rand_int(int n)
// {
// 	return( genrand_int32() % n );
// }
//
// // [0:beta)
// double rand_tau_beta(double beta)
// {
// 	return( genrand_real2()*beta );
// }
//
// // (0:l_max)
// double rand_tau_l(double l_max)
// {
// 	return( genrand_real3()*l_max );
// }


void random_permutation(int *array, int n)
{
	for(int i=0; i<n; i++)  array[i] = i;
	for(int i=0; i<n-1; i++){
		int a = rand_int(n-i);
		int temp = array[i];
		array[i] = array[i+a];
		array[i+a] = temp;
	}
}

double sgn(double x)
{
	if(x>=0)  return(+1.0);
	else  return(-1.0);
}

double pow_m1(int i)
{
	if( i%2 )  return(-1.0);
	else  return(1.0);
}

//
// 1 / (e^x + 1) = ( 1-tanh(x/2) ) / 2
//
double distrib_fermi(double x)
{
	return( 0.5 * ( 1.0 - tanh(0.5*x) ) );
}

//
// e^y / (e^x + 1)
//
double distrib_fermi_boltz(double x, double y)
{
	if(x<0)  return( exp(y) / ( exp(x) + 1.0 ) );
	else  return( exp(y-x) / ( exp(-x) + 1.0 ) );
}



int tau_order(const double *tau_array, int k, double tau)
{
	int i;

	for(i=0; i<k; i++){
		if(tau_array[i] > tau)  break;
	}

	return(i);
}
int tau_order(const std::vector<double>& array, double tau){
	return tau_order(array.data(), array.size(), tau);
}

int tau_order(const int *num_array, int k, int num)
{
	int i;

	for(i=0; i<k; i++){
		if(num_array[i] >= num)  break;
	}

	return(i);
}

int tau_position(const double *tau_array, int k, double tau)
{
	for(int i=0; i<k; i++){
		if(tau_array[i] == tau)  return(i);
	}

	printf("\n*** error in tau_position");
	exit(0);
}


//
// 0: do not update,  1: update
//
// reject if prob<0
int metropolis(double prob, unsigned long &n_accept, unsigned long &n_reject)
{
	// always accept (for test)
// 	return(1);


	if(prob<0){
		if( fabs(prob) < 1e-14 ){
			n_reject++;
			return(0);
		}

		printf("*** minus sign  (prob = %.3e)\n", prob);
// 		exit(0);
		return(0);
	}
	// (for test)
// 	if(prob<0){
// 		prob = -prob;
// 	}

	//
	// Metropolis
	//
	if( prob > 1 ){
		n_accept++;
		return(1);
	}
	else{
		// genrand_real1() generates a random number on [0,1]-real-interval
		if( prob > genrand_real1() ){
			n_accept++;
			return(1);
		}
		else{
			n_reject++;
			return(0);
		}
	}



	// Heat bath
// 	if( prob / (1.0 + prob) > genrand_real1() ){
// 		n_accept++;
// 		return(1);
// 	}
// 	else{
// 		n_reject++;
// 		return(0);
// 	}

}
//
// 0: do not update,  1: update
//
// input fabs(prob)
int metropolis_abs(double prob, unsigned long &n_accept, unsigned long &n_reject)
{
	// always accept (for test)
// 	return(1);

	//
	// Metropolis
	//
	if( prob > 1 ){
		n_accept++;
		return(1);
	}
	else{
		// genrand_real1() generates a random number on [0,1]-real-interval
		if( prob > genrand_real1() ){
			n_accept++;
			return(1);
		}
		else{
			n_reject++;
			return(0);
		}
	}
}


// return: number of data stored in params[]
int read_data(FILE *fp, double *params)
{
	const int max_str=16384;
	char str[max_str];

	do{
		if( fgets(str, max_str, fp) == NULL)  return 0;
// 		printf("%s", str);
	}while(str[0] == '#');

	int n=0;
	char *p_str = strtok(str, " ");
	while( p_str != NULL ){
		sscanf(p_str, "%lf", &params[n++]);
		p_str = strtok(NULL, " ");
	}

	return n;
}


void print_time(clock_t time1, clock_t time2)
{
	double time_dif = (double)(time2 - time1) / (double)CLOCKS_PER_SEC;

	int time_m = (int)time_dif / 60;

	double time_s = time_dif - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;

	printf("\nComputational time : ");
	if(time_h)  printf(" %dh", time_h);
	printf(" %dm %.3lfs\n", time_m, time_s);
}

void print_time(clock_t time1, clock_t time2, char *filename)
{
	double time_dif = (double)(time2 - time1) / (double)CLOCKS_PER_SEC;

	int time_m = (int)time_dif / 60;

	double time_s = time_dif - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;

	FILE *fp=fopen(filename, "a");

// 	fprintf(fp, "\nComputational time : ");
	if(time_h)  fprintf(fp, " %dh", time_h);
	fprintf(fp, " %dm %.3lfs\n", time_m, time_s);

	fclose(fp);
}

void print_time_mpi(double time_trans)
{
	int time_m = (int)time_trans / 60;

	double time_s = time_trans - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;

	printf("Transmission time : ");
	if(time_h)  printf(" %dh", time_h);
	printf(" %dm %.3lfs\n", time_m, time_s);
}

void print_time_mpi(clock_t time1, clock_t time2, double time_trans, char *filename)
{
	double time_dif = (double)(time2 - time1) / (double)CLOCKS_PER_SEC;

	int time_m = (int)time_dif / 60;

	double time_s = time_dif - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;

	FILE *fp=fopen(filename, "a");

// 	fprintf(fp, "\nComputational time : ");
	if(time_h)  fprintf(fp, " %dh", time_h);
	fprintf(fp, " %dm %.3lfs", time_m, time_s);


	int trans_m = (int)time_trans / 60;

	double trans_s = time_trans - trans_m*60;

	int trans_h = trans_m / 60;
	trans_m = trans_m % 60;

	fprintf(fp, "  (");
	if(time_h)  fprintf(fp, " %dh", trans_h);
	fprintf(fp, " %dm %.3lfs)\n", trans_m, trans_s);


	fclose(fp);
}

void sprint_time(char *str, clock_t time1)
{
	double time_dif = (double)time1 / (double)CLOCKS_PER_SEC;

	int time_m = (int)time_dif / 60;
	double time_s = time_dif - time_m*60;
	int time_h = time_m / 60;
	time_m = time_m % 60;


	char str2[100];

	sprintf(str, "");

	if(time_h){
		sprintf(str2, " %dh", time_h);
		strcat(str, str2);
	}
	sprintf(str2, " %dm %.3lfs\n", time_m, time_s);
	strcat(str, str2);
}

void sprint_time_mpi(char *str, clock_t time1, double time_trans)
{
	sprint_time(str, time1);
	str[strlen(str)-1] = '\0';  // delete '\n'

	char str2[100];

	int time_m = (int)time_trans / 60;

	double time_s = time_trans - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;

	strcat(str, "  (");
	if(time_h){
		sprintf(str2, " %dh", time_h);
		strcat(str, str2);
	}
	sprintf(str2, " %dm %.3lfs)\n", time_m, time_s);
	strcat(str, str2);
}
/*
void sprint_time_mpi(char *str, clock_t time1, double time_trans)
{
	double time_dif = (double)time1 / (double)CLOCKS_PER_SEC;

	int time_m = (int)time_dif / 60;

	double time_s = time_dif - time_m*60;

	int time_h = time_m / 60;
	time_m = time_m % 60;


	char str2[100];

	sprintf(str, "");

	if(time_h){
		sprintf(str2, " %dh", time_h);
		strcat(str, str2);
	}
	sprintf(str2, " %dm %.3lfs", time_m, time_s);
	strcat(str, str2);


	int trans_m = (int)time_trans / 60;

	double trans_s = time_trans - trans_m*60;

	int trans_h = trans_m / 60;
	trans_m = trans_m % 60;

	strcat(str, "  (");
	if(time_h){
		sprintf(str2, " %dh", trans_h);
		strcat(str, str2);
	}
	sprintf(str2, " %dm %.3lfs)\n", trans_m, trans_s);
	strcat(str, str2);
}
*/

void sprint_time_mpi(char *str, clock_t time1, clock_t time2, double time_trans)
{
	sprint_time_mpi(str, time2 - time1, time_trans);
}


void mesh_linear(double *y, double x_min, double x_max, int N)
{
	for(int i=0; i<N; i++){
		y[i] = x_min + (x_max-x_min) / (double)(N-1) * (double)i;
	}
}
void mesh_log_linear(double *y, double x_min, double x_max, int N)
{
	double x_log[N];
	double x_log_min = log(x_min), x_log_max = log(x_max);

	for(int i=0; i<N; i++){
		x_log[i] = (double)i * (x_log_max-x_log_min)/(double)(N-1) + x_log_min;
		y[i] = exp(x_log[i]);
	}
}
// inpt log10(x)
static void mesh_log_linear_2(double *y, double x_log_min, double x_log_max, int N)
{
	double x_log[N];

	for(int i=0; i<N; i++){
		x_log[i] = (double)i * (x_log_max-x_log_min)/(double)(N-1) + x_log_min;
		y[i] = pow(10, x_log[i]);
	}
}

void mesh_init(int i_mesh, double *y, double x_min, double x_max, int N)
{
	mesh_init(i_mesh, y, x_min, x_max, N, 0);
}
void mesh_init(int i_mesh, double *y, double x_min, double x_max, int N, int my_rank)
{
	void (* func_mesh[3])(double *, double, double, int)={
		mesh_linear, mesh_log_linear, mesh_log_linear_2
	};

	char str_mesh[3][128]={
		"Linear mesh",
		"Logarithmic mesh",
		"Logarithmic mesh (input log10(x))",
	};

	if(my_rank==0){
		printf("\n%s\n", str_mesh[i_mesh]);
		printf("  N=%d  MIN=%lf  MAX=%lf\n", N, x_min, x_max);
	}

	func_mesh[i_mesh](y, x_min, x_max, N);

	if(my_rank==0){
		for(int i=0; i<N; i++){
			printf(" %3d  %lf = 10 ^ %+.3lf\n", i, y[i], log10(y[i]));
		}
	}
}

// void mesh(double *y, double x_min, double x_max, int N, int i_mesh)
// {
//
// 	switch(i_mesh){
// 		case 0:
// 			mesh_linear(y, x_min, x_max, N);
// 			break;
// 		case 1:
// 			mesh_log_linear(y, x_min, x_max, N);
// 			break;
// 	}
//
// }


//=========================================================================================================
//
// non-linear mesh
// (n_tau1+1) points in the range [0:beta/2]
// shortest interval: beta/2/n_tau2
//
void tau_mesh_nonlinear_boson(double *tau, int n_tau1, int n_tau2, double beta)
{
	double delta_tau2 = beta * 0.5 / (double)n_tau2;

	double tau2[n_tau2+1];
	for(int i=0; i<=n_tau2; i++)  tau2[i] = delta_tau2 * (double)i;

	tau[0] = tau2[0];

	int R_SP = n_tau2/n_tau1;
	int i1=0, i2=0;
	for(int j=0; j<R_SP; j++){
		for(int i=0; i<n_tau1/2/R_SP; i++){
			i1 ++;
			i2 += j+1;
			tau[i1] = tau2[i2];
		}
	}
	for(int j=0; j<R_SP; j++){
		for(int i=0; i<n_tau1/2/R_SP; i++){
			i1 ++;
			i2 += j+R_SP;
			tau[i1] = tau2[i2];
		}
	}
}

// (n_tau1+1) points in the range [0:beta]
// shortest interval: beta/n_tau2
void tau_mesh_nonlinear_fermion(double *tau, int n_tau1, int n_tau2, double beta)
{
	tau_mesh_nonlinear_boson(tau, n_tau1/2, n_tau2/2, beta);

	for(int i=0; i<n_tau1/2; i++){
		tau[n_tau1-i] = beta - tau[i];
	}
}
