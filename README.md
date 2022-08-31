# cthyb

Continuous-time quantum Monte Carlo method for an impurity Anderson model. The segment algorithm in the hybridization expansion is implemented.

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

## DCore

This code can be used as the impurity solver of DCore.
A minimal input parameters are
```
[impurity_solver]
name = JO/cthyb-seg
exec_path{str} = hybqmc
MC.n_msr{int} = 10000
```
Typical choices for ``n_msr`` are

- 10000 : No enough accuracy; But, can be used to get rough convergence of the bath function. 
- 100000 : Reasonable accuracy.
- 1000000 : High accuracy; For publication.

Other parameters can be specified as follows:
In order to specify, e.g., ``rand_seed`` parameter in ``[control]`` section, add  ``control.rand_seed{int} = 0``.
