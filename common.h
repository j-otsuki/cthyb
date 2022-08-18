#ifndef _COMMON_H
#define _COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>


// MAX, MIN

template<class T>
inline T common_max(const T& a, const T& b)
{
	return ( a > b ) ? ( a ) : ( b );
}
template<class T>
inline T common_min(const T& a, const T& b)
{
	return ( a < b ) ? ( a ) : ( b );
}


// MEMORY ALLOCATION

// How to use:
// 2d
//	double **p;  // any type
//	array2_alloc(p, n1, n2);  // p[n1][n2]
//	array2_free(p, n1);
// 3d
//	double ***p;  // any type
//	array3_alloc(p, n1, n2, n3);  // p[n1][n2][n3]
//	array3_free(p, n1, n2);



template <class T>
inline void array1_alloc(T* (&p), int n1)
{
	try{
		p = new T [n1];
	}
	catch( std::bad_alloc ){
		printf("*** Failed to allocate memory\n");
		exit(0);
	}
}
template <class T>
inline void array1_free(T* p)
{
	delete [] p;
}
template <class T>
inline void array2_alloc(T** (&p), int n1, int n2)
{
	p = new T* [n1];
	for(int i=0; i<n1; i++)  array1_alloc(p[i], n2);
}
template <class T>
inline void array2_free(T** p, int n1)
{
	for(int i=0; i<n1; i++)  array1_free(p[i]);
	delete [] p;
}
template <class T>
inline void array3_alloc(T*** (&p), int n1, int n2, int n3)
{
	p = new T** [n1];
	for(int i=0; i<n1; i++)  array2_alloc(p[i], n2, n3);
}
template <class T>
inline void array3_free(T*** p, int n1, int n2)
{
	for(int i=0; i<n1; i++)  array2_free(p[i], n2);
	delete [] p;
}
template <class T>
inline void array4_alloc(T**** (&p), int n1, int n2, int n3, int n4)
{
	p = new T*** [n1];
	for(int i=0; i<n1; i++)  array3_alloc(p[i], n2, n3, n4);
}
template <class T>
inline void array4_free(T**** p, int n1, int n2, int n3)
{
	for(int i=0; i<n1; i++)  array3_free(p[i], n2, n3);
	delete [] p;
}
template <class T>
inline void array5_alloc(T***** (&p), int n1, int n2, int n3, int n4, int n5)
{
	p = new T**** [n1];
	for(int i=0; i<n1; i++)  array4_alloc(p[i], n2, n3, n4, n5);
}
template <class T>
inline void array5_free(T***** p, int n1, int n2, int n3, int n4)
{
	for(int i=0; i<n1; i++)  array4_free(p[i], n2, n3, n4);
	delete [] p;
}

template <class T>
inline void array1_zero(T* p, int n1)
{
	for(int i=0; i<n1; i++)  p[i] = 0;
}
template <class T>
inline void array2_zero(T** p, int n1, int n2)
{
	for(int i=0; i<n1; i++)  array1_zero(p[i], n2);
}
template <class T>
inline void array3_zero(T*** p, int n1, int n2, int n3)
{
	for(int i=0; i<n1; i++)  array2_zero(p[i], n2, n3);
}
template <class T>
inline void array4_zero(T**** p, int n1, int n2, int n3, int n4)
{
	for(int i=0; i<n1; i++)  array3_zero(p[i], n2, n3, n4);
}
template <class T>
inline void array5_zero(T***** p, int n1, int n2, int n3, int n4, int n5)
{
	for(int i=0; i<n1; i++)  array4_zero(p[i], n2, n3, n4, n5);
}

template <class T>
inline void array1_copy(T* p_from, T* p_to, int n1)
{
	for(int i=0; i<n1; i++)  p_to[i] = p_from[i];
}
template <class T>
inline void array2_copy(T** p_from, T** p_to, int n1, int n2)
{
	for(int i=0; i<n1; i++)  array1_copy(p_from[i], p_to[i], n2);
}
template <class T>
inline void array3_copy(T*** p_from, T*** p_to, int n1, int n2, int n3)
{
	for(int i=0; i<n1; i++)  array2_copy(p_from[i], p_to[i], n2, n3);
}
template <class T>
inline void array4_copy(T**** p_from, T**** p_to, int n1, int n2, int n3, int n4)
{
	for(int i=0; i<n1; i++)  array3_copy(p_from[i], p_to[i], n2, n3, n4);
}
template <class T>
inline void array5_copy(T***** p_from, T***** p_to, int n1, int n2, int n3, int n4, int n5)
{
	for(int i=0; i<n1; i++)  array4_copy(p_from[i], p_to[i], n2, n3, n4, n5);
}




// TIME

inline void clock2str(char *str, clock_t time_elapsed)
{
	double d_time = (double)time_elapsed / (double)CLOCKS_PER_SEC;
	
	int time_m = (int)d_time / 60;
	double time_s = d_time - time_m*60;
	int time_h = time_m / 60;
	time_m = time_m % 60;
	
// 	sprintf(str, "");
	
	if(time_h)  sprintf(str, "%dh %dm %.3lfs", time_h, time_m, time_s);
	else        sprintf(str, "%dm %.3lfs", time_m, time_s);
}


// MESH

// linear mesh
// void mesh_linear(double *x, double x_ini, double x_fin, int N);

// logarithmic mesh
// void mesh_log(double *x, double x_ini, double x_fin, int N);

// void (*const func_mesh[2])(double *, double, double, int) = {
// 	mesh_linear, mesh_log
// };


// FILE

// Read one line in the file 'fp' and store the data in 'param[]'.
// An empty line and a line starting with '#' are skipped.
// return: number of data stored in params[]
// int read_data(FILE *fp, double *params);


#endif // _COMMON_H
