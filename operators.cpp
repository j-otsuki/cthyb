/*

Continuous-time quantum Monte Carlo

written by Junya Otsuki, 2009
Dept. of Physics, Tohoku University, Sendai, Japan

*/

#include "operators.h"
#include "ct_qmc_share.h"


Operators::Operators(int max_k)
	: max_k(max_k)
	, k(0)
	, D(max_k)
	, wind(false)
	, tau1(max_k)
	, tau2(max_k)
	, beta(-1e100)
{}


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


bool Operators::is_occupied(double tau) const
{
	int i_tau1 = tau_order(tau1, tau);
	int i_tau2 = tau_order(tau2, tau);

	bool occupied;
	if( !wind ){  // no wind
		occupied = (i_tau1 != i_tau2);
	}
	else{  // wind
		occupied = (i_tau1 == i_tau2);
	}
	return occupied;
}

// total segment length; This function works for all values of k including k=0.
double Operators::length() const
{
	double l=0;
	for(int i=0; i<k; i++){
		l += tau1[i] - tau2[i];
	}
	if(wind)  l += beta;

	return(l);
}

double Operators::overlap(double tau_from, double tau_to) const
{
	if( k == 0 ){
		double l_over = 0;
		if( wind ){  // |1>
			l_over = tau_to - tau_from;
			if( l_over < 0 )  l_over += beta;
		}
		return l_over;
	}

	// double *tau1, *tau2;
	const std::vector<double> *tau1, *tau2;
	double *l1, *l2, l_comp, l_over, l_total;
	int n=0;

	if( !wind ){
		tau1 = &this->tau1;
		tau2 = &this->tau2;
		l1 = &l_over;
		l2 = &l_comp;
	}
	else{
		tau1 = &this->tau2;
		tau2 = &this->tau1;
		l1 = &l_comp;
		l2 = &l_over;
	}

	*l1 = 0;

	int i_from = tau_order(*tau2, tau_from);
	int i_to = tau_order(*tau2, tau_to);

	#if DEBUG
	printf("\ntau_overlap\n");
	printf("  i_from=%d, i_to=%d\n", i_from, i_to);
	#endif

	if(tau_from < tau_to){

		for(int i=i_from; i<i_to; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if DEBUG
			printf(" (*tau1)[%d]-tau2[%d]=%.5lf\n", i, i, (*tau1)[i] - tau2[i]);
			#endif
		}

		if(i_from){  // !=0
			double under_estim = (*tau1)[i_from-1] - tau_from;
			if( under_estim > 0 ){
				*l1 += under_estim;
				n += 1;

				#if DEBUG
				printf(" under_estim=%.5lf\n", under_estim);
				#endif
			}
		}

		if(i_to){  // !=0
			double over_estim = (*tau1)[i_to-1] - tau_to;
			if( over_estim > 0 ){
				*l1 -= over_estim;
				n -= 1;

				#if DEBUG
				printf(" over_estim=%.5lf\n", over_estim);
				#endif
			}
		}

		l_total = tau_to - tau_from;
		*l2 = l_total - *l1;
	}
	else{

		for(int i=i_from; i<k; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if DEBUG
			printf(" tau1[%d]-tau2[%d]=%.5lf\n", i, i, (*tau1)[i] - (*tau2)[i]);
			#endif
		}

		if(i_from){  // !=0
			double under_estim = (*tau1)[i_from-1] - tau_from;
			if( under_estim > 0 ){
				*l1 += under_estim;
				n += 1;

				#if DEBUG
				printf(" under_estim=%.5lf\n", under_estim);
				#endif
			}
		}

		for(int i=0; i<i_to; i++){
			*l1 += (*tau1)[i] - (*tau2)[i];
			n += 2;

			#if DEBUG
			printf(" (*tau1)[%d]-(*tau2)[%d]=%.5lf\n", i, i, (*tau1)[i] - (*tau2)[i]);
			#endif
		}

		if(i_to){  // !=0
			double over_estim = (*tau1)[i_to-1] - tau_to;
			if( over_estim > 0 ){
				*l1 -= over_estim;
				n -= 1;

				#if DEBUG
				printf(" over_estim=%.5lf\n", over_estim);
				#endif
			}
		}

		l_total = beta + tau_to - tau_from;
		*l2 = l_total - *l1;
	}

	#if DEBUG
	printf("  l_overlap=%.5lf  l_total=%.5lf\n", l_over, l_total);
	printf("  n_overlap=%d\n", n);
	#endif

	return( l_over );
}
