/*!
 \file pade_c.cpp
 \brief Definition of class Pade
*/

#include "pade.h"

using namespace std;


Pade::Pade(const vector<complex<double> > &z, const vector<complex<double> > &u)
{
	init(z, u);
}

void Pade::init(const vector<complex<double> > &z, const vector<complex<double> > &u)
{
	m_n = min(z.size(), u.size());
	m_z = z;  // copy
	m_a.resize(m_n);

	// g_n(z_i)
	vector<complex<double> > g0(u);  // copy
	vector<complex<double> > g1(m_n);

	m_a[0] = u[0];
	for(int n=1; n<m_n; n++){
		// g0[i] = g_{n-1}[i]
		// g1[i] = g_{n}[i]
		for(int i=n; i<m_n; i++)  g1[i] = ( g0[n-1] / g0[i] - 1.0 ) / (z[i] - z[n-1]);
		m_a[n] = g1[n];  // g[n][n]
		g0 = g1;
	}
}

complex<double> Pade::eval(complex<double> w)
{
	complex<double> A0 = 0.0, A1 = m_a[0], A2, B0 = 1.0, B1 = 1.0, B2;

	for(int i=1; i<m_n; i++){  // (N-1)/2
		A2 = A1 + ( w - m_z[i-1] ) * m_a[i] * A0;
		B2 = B1 + ( w - m_z[i-1] ) * m_a[i] * B0;

// 		printf("%d  %e  %e\n", i, abs(A2), abs(B2));
		A1 /= B2;
		A2 /= B2;
		B1 /= B2;
		B2 /= B2;

		A0 = A1;
		A1 = A2;
		B0 = B1;
		B1 = B2;
	}

	return( A2 / B2 );
}
