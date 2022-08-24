/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "op_determinant.h"


// OpDet::OpDet(int max_k) : max_k(max_k), __n_k(0), flag(0)
OpDet::OpDet(int max_k) : max_k(max_k), _n_k(0)
{
	mat_M.resize(max_k);
	for(int i=0; i<max_k; i++){
		mat_M[i].resize(max_k);
	}

}


double OpDet::add_lambda(double *G0_tau1, double *G0_tau2, double diag)
{
	double temp_sum = 0;
	for(int i=0; i<_n_k; i++){
		for(int j=0; j<_n_k; j++){
			temp_sum += G0_tau2[i] * mat_M[i][j] * G0_tau1[j];
		}
	}

	return( diag - temp_sum );
}

void OpDet::add_mat_M(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda)
{
	double L[_n_k], R[_n_k];
	for(int i=0; i<_n_k; i++){
		L[i] = 0;
		R[i] = 0;
		for(int l=0; l<_n_k; l++){
			L[i] += mat_M[i][l] * G0_tau1[l];
			R[i] += G0_tau2[l] * mat_M[l][i];
		}
	}
// 	for(int i=0; i<_n_k; i++)  printf("L=%8.5lf  R=%8.5lf\n", L[i], R[i]);

	double lambda_inv = 1.0/lambda;
	for(int i=0; i<_n_k; i++){
		for(int j=0; j<_n_k; j++){
			mat_M[i][j] += L[i] * R[j] * lambda_inv;
		}
	}

// 	printf("\nmat_M (k-1)\n");
// 	for(int i=0; i<_n_k; i++){
// 		for(int j=0; j<_n_k; j++){
// 			printf(" %10.5lf", mat_M[i][j]);
//
// 		}
// 		printf("\n");
// 	}

	// insert new data into row [i_tau1] & column [i_tau2]
	for(int i=0; i<_n_k; i++){
		for(int j=_n_k; j>i_tau2; j--)  mat_M[i][j] = mat_M[i][j-1];
		mat_M[i][i_tau2] = -L[i] * lambda_inv;
	}

	for(int j=0; j<=_n_k; j++){
		for(int i=_n_k; i>i_tau1; i--)  mat_M[i][j] = mat_M[i-1][j];
	}
	for(int j=0; j<_n_k; j++){
		if( j<i_tau2 )  mat_M[i_tau1][j] = -R[j] * lambda_inv;
		else  mat_M[i_tau1][j+1] = -R[j] * lambda_inv;
	}
	mat_M[i_tau1][i_tau2] = lambda_inv;

	++_n_k;
}


double OpDet::remove_lambda(int i_tau1, int i_tau2)
{
	return mat_M[i_tau1][i_tau2];
}


void OpDet::remove_mat_M(int i_tau1, int i_tau2, double lambda)
{
	double temp_M_i[_n_k-1], temp_M_j[_n_k-1];
	for(int i=0; i<_n_k-1; i++){
		if( i<i_tau1 )  temp_M_i[i] = mat_M[i][i_tau2];
		else  temp_M_i[i] = mat_M[i+1][i_tau2];

		if( i<i_tau2 )  temp_M_j[i] = mat_M[i_tau1][i];
		else  temp_M_j[i] = mat_M[i_tau1][i+1];
	}

	for(int i=0; i<_n_k; i++){
		for(int j=i_tau2; j<_n_k-1; j++)  mat_M[i][j] = mat_M[i][j+1];
	}
	for(int j=0; j<_n_k-1; j++){
		for(int i=i_tau1; i<_n_k-1; i++)  mat_M[i][j] = mat_M[i+1][j];
	}

	double lambda_inv = 1.0 / lambda;

	for(int i=0; i<_n_k-1; i++){
		for(int j=0; j<_n_k-1; j++){
			mat_M[i][j] -= temp_M_i[i] * temp_M_j[j] * lambda_inv;
		}
	}

	--_n_k;
}


double OpDet::shift1_lambda(int i_tau1, double *G0_tau1)
{
	double lambda = 0;
	for(int i=0; i<_n_k; i++){
		lambda += mat_M[i_tau1][i] * G0_tau1[i];
	}

	return(lambda);
}

void OpDet::shift1_mat_M(int i_tau1, double *G0_tau1, double lambda)
{
	double L[_n_k];
	for(int i=0; i<_n_k; i++){
		L[i] = 0;
		for(int l=0; l<_n_k; l++)  L[i] += mat_M[i][l] * G0_tau1[l];
	}

	double lambda_inv = 1.0 / lambda;

	for(int j=0; j<_n_k; j++)  mat_M[i_tau1][j] *= lambda_inv;

	for(int i=0; i<i_tau1; i++){
		for(int j=0; j<_n_k; j++)  mat_M[i][j] -= L[i] * mat_M[i_tau1][j];
	}
	for(int i=i_tau1+1; i<_n_k; i++){
		for(int j=0; j<_n_k; j++)  mat_M[i][j] -= L[i] * mat_M[i_tau1][j];
	}
}


double OpDet::shift2_lambda(int i_tau2, double *G0_tau2)
{
	double lambda = 0;
	for(int i=0; i<_n_k; i++){
		lambda += G0_tau2[i] * mat_M[i][i_tau2];
	}

	return(lambda);
}

void OpDet::shift2_mat_M(int i_tau2, double *G0_tau2, double lambda)
{
	double R[_n_k];
	for(int i=0; i<_n_k; i++){
		R[i] = 0;
		for(int l=0; l<_n_k; l++)  R[i] += G0_tau2[l] * mat_M[l][i];
	}

	double lambda_inv = 1.0 / lambda;

	for(int i=0; i<_n_k; i++)  mat_M[i][i_tau2] *= lambda_inv;

	for(int i=0; i<_n_k; i++){
		for(int j=0; j<i_tau2; j++)  mat_M[i][j] -= mat_M[i][i_tau2] * R[j];
	}
	for(int i=0; i<_n_k; i++){
		for(int j=i_tau2+1; j<_n_k; j++)  mat_M[i][j] -= mat_M[i][i_tau2] * R[j];
	}
}


double OpDet::shift12_lambda(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double diag)
{
	double L_diag = 0;
	for(int i=0; i<_n_k; i++)  L_diag += mat_M[i_tau1][i] * G0_tau1[i];

	double R_diag = 0;
	for(int j=0; j<_n_k; j++){
		if( j != i_tau1 ){
			R_diag += G0_tau2[j] * mat_M[j][i_tau2];
		}
	}


	double temp_sum = 0;
	for(int i=0; i<_n_k; i++){
		for(int j=0; j<_n_k; j++){
			if( j != i_tau1 ){
				temp_sum += G0_tau2[j] * mat_M[j][i] * G0_tau1[i];
			}
		}
	}

	double lambda = L_diag * R_diag - ( temp_sum - diag ) * mat_M[i_tau1][i_tau2];

	return(lambda);
}


void OpDet::shift12_mat_M(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda)
{
	double L[_n_k];
	for(int j=0; j<_n_k; j++){
		L[j] = 0;
		for(int i=0; i<_n_k; i++)  L[j] += mat_M[j][i] * G0_tau1[i];
	}

	double R[_n_k];
	for(int i=0; i<_n_k; i++){
		R[i] = 0;
		for(int j=0; j<_n_k; j++){
			if( j != i_tau1 ){
				R[i] += G0_tau2[j] * mat_M[j][i];
			}
		}
	}

	double temp_tau12 = mat_M[i_tau1][i_tau2];

	double temp_tau1[_n_k], temp_tau2[_n_k];
	for(int i=0; i<_n_k; i++){
		temp_tau1[i] = mat_M[i_tau1][i];
		temp_tau2[i] = mat_M[i][i_tau2];
	}

	double lambda_inv = 1.0 / lambda;

	mat_M[i_tau1][i_tau2] *= lambda_inv;

	for(int j=0; j<_n_k; j++){
		if( j != i_tau1 ){
			mat_M[j][i_tau2] = lambda_inv
			 * (L[i_tau1] * mat_M[j][i_tau2] - L[j] * temp_tau12);
		}
	}

	for(int i=0; i<_n_k; i++){
		if( i != i_tau2 ){
			mat_M[i_tau1][i] = lambda_inv
			 * (mat_M[i_tau1][i] * R[i_tau2] - temp_tau12 * R[i]);
		}
	}

	for(int j=0; j<_n_k; j++){
		if( j != i_tau1 ){
			for(int i=0; i<_n_k; i++){
				if( i != i_tau2 ){
					mat_M[j][i]
					 += mat_M[j][i_tau2] * mat_M[i_tau1][i] / mat_M[i_tau1][i_tau2]
					 - temp_tau2[j] * temp_tau1[i] / temp_tau12;
				}
			}
		}
	}
}


void OpDet::rotate_columns_upward()
{
	for(int j=0; j<_n_k; j++){
		double temp = mat_M[j][_n_k-1];

		for(int i=_n_k-1; i>0; i--)  mat_M[j][i] = mat_M[j][i-1];

		mat_M[j][0] = temp;
	}
}

void OpDet::rotate_columns_downward()
{
	for(int j=0; j<_n_k; j++){
		double temp = mat_M[j][0];

		for(int i=0; i<_n_k-1; i++)  mat_M[j][i] = mat_M[j][i+1];

		mat_M[j][_n_k-1] = temp;
	}
}

void OpDet::rotate_rows_upward()
{
	for(int j=0; j<_n_k; j++){
		double temp = mat_M[_n_k-1][j];

		for(int i=_n_k-1; i>0; i--)  mat_M[i][j] = mat_M[i-1][j];

		mat_M[0][j] = temp;
	}
}

void OpDet::rotate_rows_downward()
{
	for(int j=0; j<_n_k; j++){
		double temp = mat_M[0][j];

		for(int i=0; i<_n_k-1; i++)  mat_M[i][j] = mat_M[i+1][j];

		mat_M[_n_k-1][j] = temp;
	}
}
