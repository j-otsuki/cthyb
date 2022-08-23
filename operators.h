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
    Operators() {};
    Operators(int max_k);

    std::vector<double> tau1;  // for f-annihilation (c-creation) operator
	std::vector<double> tau2;  // for f-creation (c-annihilation) operator
    bool wind;  // true: tau1[i] < tau2[i], false: tau1[i] > tau2[i]

    OpDet D;

    int k;

    void insert_tau1(double tau);
    void insert_tau2(double tau);

    void remove_tau1(int);
    void remove_tau2(int);

    bool is_occupied(double tau) const;

    // beta must be set before length() or overlap() are called.
    void set_beta(double beta){ this->beta=beta; }
    double length() const;
    double overlap(double tau_from, double tau_to) const;

    //rotate
    void rotate_upward_tau1();
    void rotate_upward_tau2();
    void rotate_downward_tau1();
    void rotate_downward_tau2();

private:
    // int _n_k;  // number of segments
    int max_k;
    double beta;
};

#endif // _OPERATORS_H
