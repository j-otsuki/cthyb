#ifndef _CPPLAPACK_SUPPL_H
#define _CPPLAPACK_SUPPL_H

#include "cpplapack.h"


// DIAGONALIZATION ============================================================
// mat is destroyed

// dsymatrix
//! Diagonalize dsymatrix. ut * mat * u = diag(w), where ut is transpose of u
int inline CPPL_dsyev(CPPL::dsymatrix &mat, std::vector<double> &w, CPPL::dgematrix &u, CPPL::dgematrix &ut)
{
	// mat.m == mat.n
	int n = mat.n;
	std::vector< CPPL::dcovector > temp_u;
	int info = (int)mat.dsyev(w, temp_u);
	
	u.resize(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)  u(i,j) = temp_u[j](i);
	}
	ut = CPPL::t(u);
	
	return info;
}

// dgematrix
//! Diagonalize dgematrix. vr_inv * mat * vr = diag(w), where vr_inv is inverse of vr
int inline CPPL_dgeev(CPPL::dgematrix &mat, std::vector< std::complex<double> > &w, CPPL::zgematrix &vr, CPPL::zgematrix &vr_inv)
{
	// mat.m == mat.n
	int n = mat.n;
	std::vector<double> wr, wi;
	std::vector< CPPL::dcovector > temp_vrr, temp_vri;
	int info = (int)mat.dgeev(wr, wi, temp_vrr, temp_vri);
	
	w.resize(n);
	for(int i=0; i<n; i++)  w[i] = std::complex<double>(wr[i], wi[i]);
	
	vr.resize(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			vr(i,j) = std::complex<double>(temp_vrr[j](i), temp_vri[j](i));
	}
	vr_inv = CPPL::i(vr);
	
	return info;
}

// zhematrix
//! Diagonalize zhematrix. uh * mat * u = diag(w), where uh is Hermitian conjugate of u
int inline CPPL_zheev(CPPL::zhematrix &mat, std::vector<double> &w, CPPL::zgematrix &u, CPPL::zgematrix &uh)
{
	// mat.m == mat.n
	int n = mat.n;
	std::vector< CPPL::zcovector > temp_u;
	int info = (int)mat.zheev(w, temp_u);
	
	u.resize(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)  u(i,j) = temp_u[j](i);
	}
	uh = CPPL::conjt(u);
	
	return info;
}

// zgematrix
//! Diagonalize zgematrix. vr_inv * mat * vr = diag(w), where vr_inv is inverse of vr
int inline CPPL_zgeev(CPPL::zgematrix &mat, std::vector< std::complex<double> > &w, CPPL::zgematrix &vr, CPPL::zgematrix &vr_inv)
{
	// mat.m == mat.n
	int n = mat.n;
	std::vector< CPPL::zcovector > temp_vr;
	int info = (int)mat.zgeev(w, temp_vr);
	
	vr.resize(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)  vr(i,j) = temp_vr[j](i);
	}
	vr_inv = CPPL::i(vr);
	
	return info;
}

// ============================================================================
// product of a matrix and a diagonal matrix


//! dgematrix = dgematrix * diag(real)
inline CPPL::dgematrix operator*(const CPPL::dgematrix &mat, const std::vector<double> &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	CPPL::dgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = mat(i,j) * diag[j];
	}
	return r;
}
//! zgematrix = zgematrix * diag(real)
inline CPPL::zgematrix operator*(const CPPL::zgematrix &mat, const std::vector<double> &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	CPPL::zgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = mat(i,j) * diag[j];
	}
	return r;
}
//! zgematrix = zgematrix * diag(complex)
inline CPPL::zgematrix operator*(const CPPL::zgematrix &mat, const std::vector< std::complex<double> > &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	CPPL::zgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = mat(i,j) * diag[j];
	}
	return r;
}


//! dgematrix = diag(real) * dgematrix
inline CPPL::dgematrix operator*(const std::vector<double> &diag, const CPPL::dgematrix &mat)
{
#ifdef  CPPL_DEBUG
	if( mat.m != diag.size() )  exit(1);
#endif

	CPPL::dgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = diag[i] * mat(i,j);
	}
	return r;
}
//! zgematrix = diag(real) * zgematrix
inline CPPL::zgematrix operator*(const std::vector<double> &diag, const CPPL::zgematrix &mat)
{
#ifdef  CPPL_DEBUG
	if( mat.m != diag.size() )  exit(1);
#endif

	CPPL::zgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = diag[i] * mat(i,j);
	}
	return r;
}
//! zgematrix = diag(complex) * zgematrix
inline CPPL::zgematrix operator*(const std::vector< std::complex<double> > &diag, const CPPL::zgematrix &mat)
{
#ifdef  CPPL_DEBUG
	if( mat.m != diag.size() )  exit(1);
#endif

	CPPL::zgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  r(i,j) = diag[i] * mat(i,j);
	}
	return r;
}


//! dgematrix *= diag(real)
inline CPPL::dgematrix& operator*=(CPPL::dgematrix &mat, const std::vector<double> &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  mat(i,j) *= diag[j];
	}
	return mat;
}
//! zgematrix *= diag(real)
inline CPPL::zgematrix& operator*=(CPPL::zgematrix &mat, const std::vector<double> &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  mat(i,j) *= diag[j];
	}
	return mat;
}
//! zgematrix *= diag(complex)
inline CPPL::zgematrix& operator*=(CPPL::zgematrix &mat, const std::vector< std::complex<double> > &diag)
{
#ifdef  CPPL_DEBUG
	if( mat.n != diag.size() )  exit(1);
#endif

	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++)  mat(i,j) *= diag[j];
	}
	return mat;
}

// ============================================================================


/*
static double chop(double x)
{
	if( fabs(x) < 1e-14 )  x = 0;
	return x;
}

static std::complex<double> chop(std::complex<double> z)
{
	double zr = 0, zi = 0;
	if( fabs(real(z)) > 1e-14 )  zr = real(z);
	if( fabs(imag(z)) > 1e-14 )  zi = imag(z);
	return std::complex<double>(zr,zi);
}

CPPL::dgematrix CPPL_chop(CPPL::dgematrix &mat)
{
	CPPL::dgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++){
			r(i,j) =  chop( mat(i,j) );
		}
	}
	return r;
}

CPPL::zgematrix CPPL_chop(CPPL::zgematrix &mat)
{
	CPPL::zgematrix r(mat.m, mat.n);
	for(int i=0; i<mat.m; i++){
		for(int j=0; j<mat.n; j++){
			r(i,j) =  chop( mat(i,j) );
		}
	}
	return r;
}

// void CPPL_print(CPPL::zgematrix &mat)
// {
// 	for(int i=0; i<mat.m; i++){
// 		for(int j=0; j<mat.n; j++){
// 			std::cout << " " << chop( mat(i,j) );
// 		}
// 		std::cout << std::endl;
// 	}
// }
// void CPPL_print(CPPL::dgematrix &mat)
// {
// 	for(int i=0; i<mat.m; i++){
// 		for(int j=0; j<mat.n; j++){
// 			std::cout << " " << chop( mat(i,j) );
// 		}
// 		std::cout << std::endl;
// 	}
// }

*/

#endif // _CPPLAPACK_SUPPL_H
