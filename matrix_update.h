/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _MATRIX_UPDATE_H
#define _MATRIX_UPDATE_H

#include "ct_qmc_share.h"
//
// struct cond_op
// N_K
// LAPACK
//


//
// G0_tau1[] = G0( tau2[], tau1' )
// G0_tau2[] = G0( tau2', tau1[] )
// diag = G0( tau2', tau1' )
//
double add_lambda(struct cond_op &, double *G0_tau1, double *G0_tau2, double diag);
void add_mat_M(struct cond_op &, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda);

double add2_lambda(struct cond_op &F, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2]);
//
// add two rows and two columns at the end of the matrix M
void add2_mat_M(struct cond_op &F, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2]);
//
// insert two rows and two columns at i_tau
void add2_mat_M(struct cond_op &F, int i_tau, double G0_tau1[2][N_K], double G0_tau2[2][N_K], double lambda[2][2]);

double add2_lambda(struct cond_op &F, double **G0_tau1, double **G0_tau2, double **lambda, int n, long *ipiv);
void add2_mat_M(struct cond_op &F, double **G0_tau1, double **G0_tau2, double **lambda, int n, long *ipiv);

void remove_mat_M(struct cond_op &, int i_tau1, int i_tau2, double);

void remove2_mat_M(struct cond_op &F, int i_tau1, int i_tau2, double lambda[2][2]);

void remove2_mat_M(struct cond_op &F, double **lambda, int n, long *ipiv);


double determinant(double **mat, int N, long *ipiv);
void inverse_matrix(double **mat, int N, long *ipiv);


// shift tau1
double shift1_lambda(struct cond_op &, int i_tau1, double *G0_tau1);
void shift1_mat_M(struct cond_op &, int i_tau1, double *G0_tau1, double lambda);

// shift tau2
double shift2_lambda(struct cond_op &, int i_tau2, double *G0_tau2);
void shift2_mat_M(struct cond_op &, int i_tau2, double *G0_tau2, double lambda);

//shift tau1 & tau2
double shift12_lambda(struct cond_op &, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double diag);
void shift12_mat_M(struct cond_op &, int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda);


//
//rotate
//
void rotate_upward(double *tau, int k);
void rotate_downward(double *tau, int k);

void rotate_upward(int *alpha, int k);
void rotate_downward(int *alpha, int k);

void rotate_columns_upward(double mat[N_K][N_K], int k);
void rotate_columns_downward(double mat[N_K][N_K], int k);

void rotate_rows_upward(double mat[N_K][N_K], int k);
void rotate_rows_downward(double mat[N_K][N_K], int k);

//
// exchange
//
void exchange_tau(double *tau, int i1, int i2);
void exchange_columns(double mat[N_K][N_K], int i1, int i2, int k);
void exchange_rows(double mat[N_K][N_K], int i1, int i2, int k);


#endif // _MATRIX_UPDATE_H
