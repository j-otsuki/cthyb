#
# Makefile
#
#==============================================================================
# The following macros may be redefined in the file 'make.inc.$(HOSTNAME)'.
# If the file does not exist, or if definitions are not given in the file, the
# following definitions are used by default.

CC = mpicxx
CFLAGS = -O2 -fpermissive -DMPICH_IGNORE_CXX_SEEK -I $(HOME)/local/include
LDFLAGS = -L $(HOME)/local/lib

LIB_GSL = -lgsl -lgslcblas
LIB_LAPACK = -llapack -lblas -mkl
LIB_FFTW = -lfftw3
LIB_ARPACK = -larpack -lgfortran

-include make.inc.$(HOSTNAME)
# This command is ignored if the file does not exist.

#==============================================================================

.cpp.o:
	$(CC) $(CFLAGS) $*.cpp -c

CTQMC = libctqmc.a
CTQMCOBJ = ct_qmc_share.o matrix_update.o green_func_0.o dmft.o mt.o fft.o pade.o erfc.o gtau.o

$(CTQMC): $(CTQMCOBJ)
	ar rc $(CTQMC) $(CTQMCOBJ)
	ranlib $(CTQMC)

LIB_CTQMC = $(CTQMC) $(LIB_GSL) $(LIB_LAPACK)
LIB4CTQMC = $(LIB_GSL) $(LIB_LAPACK)

#==============================================================================

hyb_qmc: hyb_qmc_main.o hyb_qmc.o $(CTQMC)
	$(CC) $(LDFLAGS) $^ $(LIB4CTQMC) -o $@.out

sb_qmc: sb_qmc_main.o sb_qmc.o $(CTQMC)
	$(CC) $(LDFLAGS) $^ $(LIB4CTQMC) $(LIB_FFTW) -o $@.out


#==============================================================================

clean:
	rm -f *.o *.out a.exe *.a *.stackdump
