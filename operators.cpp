/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "operators.h"



// Operators::Operators(int max_k) : max_k(max_k), _n_k(0), D(max_k), flag(1)
Operators::Operators(int max_k) : max_k(max_k), k(0), D(max_k), flag(1)
{
	tau1.resize(max_k);
	tau2.resize(max_k);

}


static void rotate_upward(std::vector<double>& tau, int k)
{
	double temp = tau[k-1];

	for(int i=k-1; i>0; i--)  tau[i] = tau[i-1];

	tau[0] = temp;
}

static void rotate_downward(std::vector<double>& tau, int k)
{
	double temp = tau[0];

	for(int i=0; i<k-1; i++)  tau[i] = tau[i+1];

	tau[k-1] = temp;
}

void Operators::rotate_upward_tau1()
{
	rotate_upward(tau1, k);
}
void Operators::rotate_upward_tau2()
{
	rotate_upward(tau2, k);
}
void Operators::rotate_downward_tau1()
{
	rotate_downward(tau1, k);
}
void Operators::rotate_downward_tau2()
{
	rotate_downward(tau2, k);
}
