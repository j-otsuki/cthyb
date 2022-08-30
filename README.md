# cthyb

## Requirement

- MPI
- GSL
- Boost (header only)
- GGTW

## Install

Clone the repository by
```
git clone git@github.com:j-otsuki/cthyb.git
```
Move into an empty directory for build:
```
mkdir cthyb.build
cd cthyb.build
```
Now, build the project using ``cmake`` command 
```
cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local ../cthyb
make
make install
```
If succedded, an executable ``hybqmc`` is installed to $HOME/local directory.

Intel compiler can be used by
```
CC=icc CXX=icpc cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local ../cthyb
```

## Sample(s)

Samples are provided in samples/ directory.
How to run: 
```
mpirun -np 2 hybqmc params.ini
```
Results can be plotted using the gnuplot script:
```
gnuplot *plt
```
