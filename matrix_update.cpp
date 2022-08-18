/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "matrix_update.h"

double add_lambda(struct cond_op &F, double *G0_tau1, double *G0_tau2, double diag)
{
	double temp_sum = 0;
	for(int i=0; i<F.k; i++){
		for(int j=0; j<F.k; j++){
			temp_sum += G0_tau2[i] * F.mat_M[i][j] * G0_tau1[j];
		}
	}
	
	return( diag - temp_sum );
}

void add_mat_M(struct cond_op &F, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda)
{
	double L[F.k], R[F.k];
	for(int i=0; i<F.k; i++){
		L[i] = 0;
		R[i] = 0;
		for(int l=0; l<F.k; l++){
			L[i] += F.mat_M[i][l] * G0_tau1[l];
			R[i] += G0_tau2[l] * F.mat_M[l][i];
		}
	}
// 	for(int i=0; i<F.k; i++)  printf("L=%8.5lf  R=%8.5lf\n", L[i], R[i]);
	
	double lambda_inv = 1.0/lambda;
	for(int i=0; i<F.k; i++){
		for(int j=0; j<F.k; j++){
			F.mat_M[i][j] += L[i] * R[j] * lambda_inv;
		}
	}
	
// 	printf("\nmat_M (k-1)\n");
// 	for(int i=0; i<F.k; i++){
// 		for(int j=0; j<F.k; j++){
// 			printf(" %10.5lf", F.mat_M[i][j]);
// 			
// 		}
// 		printf("\n");
// 	}
	
	// insert new data into row [i_tau1] & column [i_tau2]
	for(int i=0; i<F.k; i++){
		for(int j=F.k; j>i_tau2; j--)  F.mat_M[i][j] = F.mat_M[i][j-1];
		F.mat_M[i][i_tau2] = -L[i] * lambda_inv;
	}
	
	for(int j=0; j<=F.k; j++){
		for(int i=F.k; i>i_tau1; i--)  F.mat_M[i][j] = F.mat_M[i-1][j];
	}
	for(int j=0; j<F.k; j++){
		if( j<i_tau2 )  F.mat_M[i_tau1][j] = -R[j] * lambda_inv;
		else  F.mat_M[i_tau1][j+1] = -R[j] * lambda_inv;
	}
	F.mat_M[i_tau1][i_tau2] = lambda_inv;
}


double add2_lambda(struct cond_op &F, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2])
{
	int n=2;
	
	for(int q1=0; q1<n; q1++){
		for(int q2=0; q2<n; q2++){
			
			double temp_sum = 0;
			for(int i=0; i<F.k; i++){
				for(int j=0; j<F.k; j++){
					temp_sum += G0_tau2[q1][i] * F.mat_M[i][j] * G0_tau1[q2][j];
				}
			}
			
			lambda[q1][q2] -= temp_sum;
		}
	}
	
	return( lambda[0][0] * lambda[1][1] - lambda[0][1] * lambda[1][0] );
}

void add2_mat_M(struct cond_op &F, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2])
{
	int n=2;
	
	double L[n][F.k], R[n][F.k];
	for(int q=0; q<n; q++){
		for(int i=0; i<F.k; i++){
			L[q][i] = 0;
			R[q][i] = 0;
			for(int l=0; l<F.k; l++){
				L[q][i] += F.mat_M[i][l] * G0_tau1[q][l];
				R[q][i] += G0_tau2[q][l] * F.mat_M[l][i];
			}
		}
	}
	
	double det_inv = 1.0 / (lambda[0][0] * lambda[1][1] - lambda[0][1] * lambda[1][0]);
	double lambda_inv[2][2];
	lambda_inv[0][0] = lambda[1][1] * det_inv;
	lambda_inv[1][1] = lambda[0][0] * det_inv;
	lambda_inv[0][1] = -lambda[0][1] * det_inv;
	lambda_inv[1][0] = -lambda[1][0] * det_inv;
	
	for(int i=0; i<F.k; i++){
		for(int j=0; j<F.k; j++){
			for(int q1=0; q1<n; q1++){
				for(int q2=0; q2<n; q2++){
					F.mat_M[i][j] += L[q1][i] * lambda_inv[q1][q2] * R[q2][j];
				}
			}
		}
	}
	
	// add new data
	
	for(int q1=0; q1<n; q1++){
		for(int i=0; i<F.k; i++){
			double temp=0;
			
			for(int q2=0; q2<n; q2++)  temp += L[q2][i] * lambda_inv[q2][q1];
			
			F.mat_M[i][F.k+q1] = -temp;
		}
	}
	
	for(int q1=0; q1<n; q1++){
		for(int i=0; i<F.k; i++){
			double temp=0;
			
			for(int q2=0; q2<n; q2++)  temp += lambda_inv[q1][q2] * R[q2][i];
			
			F.mat_M[F.k+q1][i] = -temp;
		}
	}
	
	for(int q1=0; q1<n; q1++){
		for(int q2=0; q2<n; q2++){
			F.mat_M[F.k+q1][F.k+q2] = lambda_inv[q1][q2];
		}
	}
}

void add2_mat_M(struct cond_op &F, int i_tau, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2])
{
	
	add2_mat_M(F, G0_tau1, G0_tau2, lambda);
	
	for(int i=0; i<F.k+2; i++){
		
		double temp1 = F.mat_M[i][F.k];
		double temp2 = F.mat_M[i][F.k+1];
		
		for(int j=F.k+1; j>i_tau+1; j--){
			F.mat_M[i][j] = F.mat_M[i][j-2];
		}
		
		F.mat_M[i][i_tau] = temp1;
		F.mat_M[i][i_tau+1] = temp2;
	}
	
	
	double temp1[F.k+2], temp2[F.k+2];
	for(int j=0; j<F.k+2; j++){
		temp1[j] = F.mat_M[F.k][j];
		temp2[j] = F.mat_M[F.k+1][j];
	}
	
	for(int i=F.k+1; i>i_tau+1; i--){
		for(int j=0; j<F.k+2; j++){
			F.mat_M[i][j] = F.mat_M[i-2][j];
		}
	}
	
	for(int j=0; j<F.k+2; j++){
		F.mat_M[i_tau][j] = temp1[j];
		F.mat_M[i_tau+1][j] = temp2[j];
	}
}


//
// return determinant of matrix lambda
// matrix lambda stores L and U matrices resultant from LU factorization of matrix lambda
// array ipiv stores pivot indices in the LU factorization
//
double add2_lambda(struct cond_op &F, double **G0_tau1, double **G0_tau2, double **lambda, int n, long *ipiv)
{
	if(n==0)  return(1.0);
	if(n==1){
		lambda[0][0] = add_lambda(F, G0_tau1[0], G0_tau2[0], lambda[0][0]);
		return( lambda[0][0] );
	}
	
	for(int q1=0; q1<n; q1++){
		for(int q2=0; q2<n; q2++){
			
			double temp_sum = 0;
			for(int i=0; i<F.k; i++){
				for(int j=0; j<F.k; j++){
					temp_sum += G0_tau2[q1][i] * F.mat_M[i][j] * G0_tau1[q2][j];
				}
			}
			
			lambda[q1][q2] -= temp_sum;
		}
	}
	
	return( determinant(lambda, n, ipiv) );
}

void add2_mat_M(struct cond_op &F, double **G0_tau1, double **G0_tau2, double **lambda, int n, long *ipiv)
{
	if(n==1){
		add_mat_M(F, F.k, F.k, G0_tau1[0], G0_tau2[0], lambda[0][0]);
		return;
	}
	
	double L[n][F.k], R[n][F.k];
	for(int q=0; q<n; q++){
		for(int i=0; i<F.k; i++){
			L[q][i] = 0;
			R[q][i] = 0;
			for(int l=0; l<F.k; l++){
				L[q][i] += F.mat_M[i][l] * G0_tau1[q][l];
				R[q][i] += G0_tau2[q][l] * F.mat_M[l][i];
			}
		}
	}
	
	inverse_matrix(lambda, n, ipiv);
	
	for(int i=0; i<F.k; i++){
		for(int j=0; j<F.k; j++){
			for(int q1=0; q1<n; q1++){
				for(int q2=0; q2<n; q2++){
					F.mat_M[i][j] += L[q1][i] * lambda[q1][q2] * R[q2][j];
				}
			}
		}
	}
	
	// add new data
	
	for(int q1=0; q1<n; q1++){
		for(int i=0; i<F.k; i++){
			double temp=0;
			
			for(int q2=0; q2<n; q2++)  temp += L[q2][i] * lambda[q2][q1];
			
			F.mat_M[i][F.k+q1] = -temp;
		}
	}
	
	for(int q1=0; q1<n; q1++){
		for(int i=0; i<F.k; i++){
			double temp=0;
			
			for(int q2=0; q2<n; q2++)  temp += lambda[q1][q2] * R[q2][i];
			
			F.mat_M[F.k+q1][i] = -temp;
		}
	}
	
	for(int q1=0; q1<n; q1++){
		for(int q2=0; q2<n; q2++){
			F.mat_M[F.k+q1][F.k+q2] = lambda[q1][q2];
		}
	}
}



void remove_mat_M(struct cond_op &F, int i_tau1, int i_tau2, double lambda)
{
	double temp_M_i[F.k-1], temp_M_j[F.k-1];
	for(int i=0; i<F.k-1; i++){
		if( i<i_tau1 )  temp_M_i[i] = F.mat_M[i][i_tau2];
		else  temp_M_i[i] = F.mat_M[i+1][i_tau2];
		
		if( i<i_tau2 )  temp_M_j[i] = F.mat_M[i_tau1][i];
		else  temp_M_j[i] = F.mat_M[i_tau1][i+1];
	}
	
	for(int i=0; i<F.k; i++){
		for(int j=i_tau2; j<F.k-1; j++)  F.mat_M[i][j] = F.mat_M[i][j+1];
	}
	for(int j=0; j<F.k-1; j++){
		for(int i=i_tau1; i<F.k-1; i++)  F.mat_M[i][j] = F.mat_M[i+1][j];
	}
	
	double lambda_inv = 1.0 / lambda;
	
	for(int i=0; i<F.k-1; i++){
		for(int j=0; j<F.k-1; j++){
			F.mat_M[i][j] -= temp_M_i[i] * temp_M_j[j] * lambda_inv;
		}
	}
}


void remove2_mat_M(struct cond_op &F, int i_tau1, int i_tau2, double lambda[2][2])
{
	int n=2;
	
	double det_inv = 1.0 / (lambda[0][0] * lambda[1][1] - lambda[0][1] * lambda[1][0]);
	double lambda_inv[2][2];
	lambda_inv[0][0] = lambda[1][1] * det_inv;
	lambda_inv[1][1] = lambda[0][0] * det_inv;
	lambda_inv[0][1] = -lambda[0][1] * det_inv;
	lambda_inv[1][0] = -lambda[1][0] * det_inv;
	
	
	double temp_M_tau1[F.k], temp_M_tau2[F.k];
	for(int j=0; j<F.k; j++){
		temp_M_tau1[j] = F.mat_M[i_tau1][j];
		temp_M_tau2[j] = F.mat_M[i_tau2][j];
	}
	
	for(int i=i_tau1; i<i_tau2-1; i++){
		for(int j=0; j<F.k; j++){
			F.mat_M[i][j] = F.mat_M[i+1][j];
		}
	}
	for(int i=i_tau2-1; i<F.k-2; i++){
		for(int j=0; j<F.k; j++){
			F.mat_M[i][j] = F.mat_M[i+2][j];
		}
	}
	
	for(int j=0; j<F.k; j++){
		F.mat_M[F.k-2][j] = temp_M_tau1[j];
		F.mat_M[F.k-1][j] = temp_M_tau2[j];
	}
	
	
	for(int j=0; j<F.k; j++){
		temp_M_tau1[j] = F.mat_M[j][i_tau1];
		temp_M_tau2[j] = F.mat_M[j][i_tau2];
	}
	
	for(int j=0; j<F.k; j++){
		for(int i=i_tau1; i<i_tau2-1; i++){
			F.mat_M[j][i] = F.mat_M[j][i+1];
		}
		for(int i=i_tau2-1; i<F.k-2; i++){
			F.mat_M[j][i] = F.mat_M[j][i+2];
		}
	}
	
	for(int j=0; j<F.k; j++){
		F.mat_M[j][F.k-2] = temp_M_tau1[j];
		F.mat_M[j][F.k-1] = temp_M_tau2[j];
	}
	
	
	for(int i=0; i<F.k-n; i++){
		for(int j=0; j<F.k-n; j++){
			for(int q1=0; q1<n; q1++){
				for(int q2=0; q2<n; q2++){
					F.mat_M[i][j] -= F.mat_M[i][F.k-n+q1] * lambda_inv[q1][q2] * F.mat_M[F.k-n+q2][j];
				}
			}
		}
	}
}

void remove2_mat_M(struct cond_op &F, double **lambda, int n, long *ipiv)
{
	inverse_matrix(lambda, n, ipiv);
	
	for(int i=0; i<F.k-n; i++){
		for(int j=0; j<F.k-n; j++){
			for(int q1=0; q1<n; q1++){
				for(int q2=0; q2<n; q2++){
					F.mat_M[i][j] -= F.mat_M[i][F.k-n+q1] * lambda[q1][q2] * F.mat_M[F.k-n+q2][j];
				}
			}
		}
	}
}


//
// return determinant of matrix mat
// matrix mat stores L and U matrices resultant from LU factorization of original matrix mat
//
double determinant(double **mat, int N, long *ipiv)
{
	#ifdef _LAPACK
	
	long int n=N, lda=N, info, lwork=N;
	double A[N*N], work[N];
	
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			A[i+j*N] = mat[i][j];
		}
	}
	
	dgetrf_(&n, &n, A, &lda, ipiv, &info);
	
	double det=1.0;
	for(int i=0; i<N-1; i++){
		if(ipiv[i] != i+1)  det = -det;
	}
	
	for(int i=0; i<N; i++)  det *= A[i+N*i];
	
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			mat[i][j] = A[i+j*N];
		}
	}
	
	return(det);
	
	#else
	
	return(0.0);
	
	#endif // _LAPACK
}

void inverse_matrix(double **mat, int N, long *ipiv)
{
	#ifdef _LAPACK
	
	long int n=N, lda=N, info, lwork=N;
	double A[N*N], work[N];
	
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			A[i+j*N] = mat[i][j];
		}
	}
	
	dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
	
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			mat[i][j] = A[i+j*N];
		}
	}
	
	#endif // _LAPACK
}



double shift1_lambda(struct cond_op &F, int i_tau1, double *G0_tau1)
{
	double lambda = 0;
	for(int i=0; i<F.k; i++){
		lambda += F.mat_M[i_tau1][i] * G0_tau1[i];
	}
	
	return(lambda);
}

void shift1_mat_M(struct cond_op &F, int i_tau1, double *G0_tau1, double lambda)
{
	double L[F.k];
	for(int i=0; i<F.k; i++){
		L[i] = 0;
		for(int l=0; l<F.k; l++)  L[i] += F.mat_M[i][l] * G0_tau1[l];
	}
	
	double lambda_inv = 1.0 / lambda;
	
	for(int j=0; j<F.k; j++)  F.mat_M[i_tau1][j] *= lambda_inv;
	
	for(int i=0; i<i_tau1; i++){
		for(int j=0; j<F.k; j++)  F.mat_M[i][j] -= L[i] * F.mat_M[i_tau1][j];
	}
	for(int i=i_tau1+1; i<F.k; i++){
		for(int j=0; j<F.k; j++)  F.mat_M[i][j] -= L[i] * F.mat_M[i_tau1][j];
	}
}


double shift2_lambda(struct cond_op &F, int i_tau2, double *G0_tau2)
{
	double lambda = 0;
	for(int i=0; i<F.k; i++){
		lambda += G0_tau2[i] * F.mat_M[i][i_tau2];
	}
	
	return(lambda);
}

void shift2_mat_M(struct cond_op &F, int i_tau2, double *G0_tau2, double lambda)
{
	double R[F.k];
	for(int i=0; i<F.k; i++){
		R[i] = 0;
		for(int l=0; l<F.k; l++)  R[i] += G0_tau2[l] * F.mat_M[l][i];
	}
	
	double lambda_inv = 1.0 / lambda;
	
	for(int i=0; i<F.k; i++)  F.mat_M[i][i_tau2] *= lambda_inv;
	
	for(int i=0; i<F.k; i++){
		for(int j=0; j<i_tau2; j++)  F.mat_M[i][j] -= F.mat_M[i][i_tau2] * R[j];
	}
	for(int i=0; i<F.k; i++){
		for(int j=i_tau2+1; j<F.k; j++)  F.mat_M[i][j] -= F.mat_M[i][i_tau2] * R[j];
	}
}


double shift12_lambda(struct cond_op &F, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double diag)
{
	double L_diag = 0;
	for(int i=0; i<F.k; i++)  L_diag += F.mat_M[i_tau1][i] * G0_tau1[i];
	
	double R_diag = 0;
	for(int j=0; j<F.k; j++){
		if( j != i_tau1 ){
			R_diag += G0_tau2[j] * F.mat_M[j][i_tau2];
		}
	}
	
	
	double temp_sum = 0;
	for(int i=0; i<F.k; i++){
		for(int j=0; j<F.k; j++){
			if( j != i_tau1 ){
				temp_sum += G0_tau2[j] * F.mat_M[j][i] * G0_tau1[i];
			}
		}
	}
	
	double lambda = L_diag * R_diag - ( temp_sum - diag ) * F.mat_M[i_tau1][i_tau2];
	
	return(lambda);
}


void shift12_mat_M(struct cond_op &F, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda)
{
	double L[F.k];
	for(int j=0; j<F.k; j++){
		L[j] = 0;
		for(int i=0; i<F.k; i++)  L[j] += F.mat_M[j][i] * G0_tau1[i];
	}
	
	double R[F.k];
	for(int i=0; i<F.k; i++){
		R[i] = 0;
		for(int j=0; j<F.k; j++){
			if( j != i_tau1 ){
				R[i] += G0_tau2[j] * F.mat_M[j][i];
			}
		}
	}
	
	double temp_tau12 = F.mat_M[i_tau1][i_tau2];
	
	double temp_tau1[F.k], temp_tau2[F.k];
	for(int i=0; i<F.k; i++){
		temp_tau1[i] = F.mat_M[i_tau1][i];
		temp_tau2[i] = F.mat_M[i][i_tau2];
	}
	
	double lambda_inv = 1.0 / lambda;
	
	F.mat_M[i_tau1][i_tau2] *= lambda_inv;
	
	for(int j=0; j<F.k; j++){
		if( j != i_tau1 ){
			F.mat_M[j][i_tau2] = lambda_inv
			 * (L[i_tau1] * F.mat_M[j][i_tau2] - L[j] * temp_tau12);
		}
	}
	
	for(int i=0; i<F.k; i++){
		if( i != i_tau2 ){
			F.mat_M[i_tau1][i] = lambda_inv
			 * (F.mat_M[i_tau1][i] * R[i_tau2] - temp_tau12 * R[i]);
		}
	}
	
	for(int j=0; j<F.k; j++){
		if( j != i_tau1 ){
			for(int i=0; i<F.k; i++){
				if( i != i_tau2 ){
					F.mat_M[j][i]
					 += F.mat_M[j][i_tau2] * F.mat_M[i_tau1][i] / F.mat_M[i_tau1][i_tau2]
					 - temp_tau2[j] * temp_tau1[i] / temp_tau12;
				}
			}
		}
	}
}


//=========================================================================================================

//
//rotate
//
void rotate_upward(double *tau, int k)
{
	double temp = tau[k-1];
	
	for(int i=k-1; i>0; i--)  tau[i] = tau[i-1];
	
	tau[0] = temp;
}
void rotate_upward(int *alpha, int k)
{
	int temp = alpha[k-1];
	
	for(int i=k-1; i>0; i--)  alpha[i] = alpha[i-1];
	
	alpha[0] = temp;
}

void rotate_downward(double *tau, int k)
{
	double temp = tau[0];
	
	for(int i=0; i<k-1; i++)  tau[i] = tau[i+1];
	
	tau[k-1] = temp;
}
void rotate_downward(int *alpha, int k)
{
	int temp = alpha[0];
	
	for(int i=0; i<k-1; i++)  alpha[i] = alpha[i+1];
	
	alpha[k-1] = temp;
}

void rotate_columns_upward(double mat[N_K][N_K], int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[j][k-1];
		
		for(int i=k-1; i>0; i--)  mat[j][i] = mat[j][i-1];
		
		mat[j][0] = temp;
	}
}

void rotate_columns_downward(double mat[N_K][N_K], int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[j][0];
		
		for(int i=0; i<k-1; i++)  mat[j][i] = mat[j][i+1];
		
		mat[j][k-1] = temp;
	}
}

void rotate_rows_upward(double mat[N_K][N_K], int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[k-1][j];
		
		for(int i=k-1; i>0; i--)  mat[i][j] = mat[i-1][j];
		
		mat[0][j] = temp;
	}
}

void rotate_rows_downward(double mat[N_K][N_K], int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[0][j];
		
		for(int i=0; i<k-1; i++)  mat[i][j] = mat[i+1][j];
		
		mat[k-1][j] = temp;
	}
}

void exchange_tau(double *tau, int i1, int i2)
{
	double temp = tau[i1];
	tau[i1] = tau[i2];
	tau[i2] = temp;
}

void exchange_columns(double mat[N_K][N_K], int i1, int i2, int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[j][i1];
		mat[j][i1] = mat[j][i2];
		mat[j][i2] = temp;
	}
}

void exchange_rows(double mat[N_K][N_K], int i1, int i2, int k)
{
	for(int j=0; j<k; j++){
		double temp = mat[i1][j];
		mat[i1][j] = mat[i2][j];
		mat[i2][j] = temp;
	}
}
