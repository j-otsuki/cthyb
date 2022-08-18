#ifndef _MT_H
#define _MT_H

#define MT_MPI

#ifdef MT_MPI
#include <mpi.h>
#endif // MT_MPI


/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);

/* initialized with time (return the seed) */
unsigned long mt_init_by_time();
unsigned long mt_init_by_time_mpi();  // MPI version
void mt_init_by_time(unsigned long *seed);  // initialized with time if *seed=0
void mt_init_by_time_mpi(unsigned long *seed);

// void mt_init_by_time(void);
// void mt_init_by_time_mpi(void);  // MPI version


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
inline long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
inline double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
inline double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
inline double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
inline double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

// /* generates a random number on [0,0x7fffffff]-interval */
// long genrand_int31(void);
// 
// /* generates a random number on [0,1]-real-interval */
// double genrand_real1(void);
// 
// /* generates a random number on [0,1)-real-interval */
// double genrand_real2(void);
// 
// /* generates a random number on (0,1)-real-interval */
// double genrand_real3(void);
// 
// /* generates a random number on [0,1) with 53-bit resolution*/
// double genrand_res53(void);


#endif // _MT_H
