# cthyb

Continuous-time quantum Monte Carlo method for an impurity Anderson model. The segment algorithm in the hybridization expansion is implemented.

## Features

To be updated.

## Requirement

- MPI
- GSL
- Boost (header only)
- FFTW

## Installation

Clone the repository by
```
$ git clone git@github.com:j-otsuki/cthyb.git
```
Move into an empty directory for build:
```
$ mkdir cthyb.build
$ cd cthyb.build
```
Now, build the project using ``cmake`` command 
```
$ cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local ../cthyb
$ make
$ make install
```
If succedded, an executable ``hybqmc`` is installed to $HOME/local/bin directory.

If one wants to use the Intel compiler, use ``CC`` and ``CXX`` environment variable as follows:
```
$ CC=icc CXX=icpc cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local ../cthyb
```

One can check if the executable is properly installed by
```
$ hybqmc --version
hybqmc version 0.2.0
```


## Samples

Samples are provided in samples/ directory.
How to run: 
```
$ mpirun -np 2 hybqmc params.ini
```
Results can be plotted using the gnuplot script:
```
$ gnuplot *plt
```

## DCore

This code can be used as an impurity solver of DCore.
A minimal input parameters for DCore are
```
[impurity_solver]
name = JO/cthyb-seg
exec_path{str} = hybqmc
MC.n_msr{int} = 10000
```
The total number of samples are ``n_msr*n_bin``, where ``n_bin=10`` in the default setting.
Typical choices for ``n_msr`` are

- 10000 : No enough accuracy; But, can be used to get rough convergence of the bath function. 
- 100000 : Reasonable accuracy.
- 1000000 : High accuracy; For publication.

Other optional parameters can be specified as well:
In order to specify, e.g., ``rand_seed`` parameter in ``[control]`` section, add  ``control.rand_seed{int} = 0``.
