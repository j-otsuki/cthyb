/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#ifndef _OPERATORS_H
#define _OPERATORS_H

#include "op_determinant.h"
#include <vector>


class Operators{
public:
    Operators(int max_k);

    std::vector<double> tau1;  // for f-annihilation (c-creation) operator
	std::vector<double> tau2;  // for f-creation (c-annihilation) operator
	int flag;  // 0: tau1[i] < tau2[i],  1: tau1[i] > tau2[i]

    OpDet D;

    int k;
    // int n_k() { return _n_k; };

    void insert_tau1(double tau);
    void insert_tau2(double tau);

    void remove_tau1(int);
    void remove_tau2(int);

    double length();

    //rotate
    void rotate_upward_tau1();
    void rotate_upward_tau2();
    void rotate_downward_tau1();
    void rotate_downward_tau2();

private:
    // int _n_k;  // number of segments
    int max_k;
};


int tau_order(std::vector<double>, double);


#endif // _OPERATORS_H
