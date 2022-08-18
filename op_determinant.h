/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _OP_DETERMINANT_H
#define _OP_DETERMINANT_H

#include <vector>


class OpDet{
public:
    OpDet() {};
    OpDet(int max_k);

    // std::vector<double> tau1;  // for f-annihilation (c-creation) operator
	// std::vector<double> tau2;  // for f-creation (c-annihilation) operator
	std::vector< std::vector<double> > mat_M;  // [tau1][tau2]
	// int flag;  // 0: tau1[i] < tau2[i],  1: tau1[i] > tau2[i]

    int n_k() { return _n_k; };

    // add a segment
    //   G0_tau1[] = G0( tau2[], tau1' )
    //   G0_tau2[] = G0( tau2', tau1[] )
    //   diag = G0( tau2', tau1' )
    double add_lambda(double *G0_tau1, double *G0_tau2, double diag);
    void add_mat_M(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda);

    // remove a segment
    double remove_lambda(int i_tau1, int i_tau2);
    void remove_mat_M(int i_tau1, int i_tau2, double);

    // shift tau1
    double shift1_lambda(int i_tau1, double *G0_tau1);
    void shift1_mat_M(int i_tau1, double *G0_tau1, double lambda);

    // shift tau2
    double shift2_lambda(int i_tau2, double *G0_tau2);
    void shift2_mat_M(int i_tau2, double *G0_tau2, double lambda);

    // shift tau1 & tau2
    double shift12_lambda(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double diag);
    void shift12_mat_M(int i_tau1, int i_tau2, double *G0_tau1, double *G0_tau2, double lambda);

    // rotate matrix
    void rotate_columns_upward();
    void rotate_columns_downward();
    void rotate_rows_upward();
    void rotate_rows_downward();

private:
	int _n_k;  // number of segments
    int max_k;
};

#endif // _OP_DETERMINANT_H
