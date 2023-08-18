# Stability-FHN-Oscillatory-Pulses
Matlab code to reproduce computations in [Carter, de Rijk, and Sandstede 2016].

# Instructions for Use

This repository provides Matlab code to reproduce the numerical results presented in Figures 6 and 7 in [P Carter, B de Rijk, and B Sandstede. Stability of traveling pulses with oscillatory tails in the FitzHugh-Nagumo system. Journal of Nonlinear Science 26 (2016) 1369–1444](http://dx.doi.org/10.1007/s00332-016-9308-7)). Details about the methods can be found in this paper.

## Matlab Scripts 

The code consists of three Matlab scripts:

`solve_FHN_w.m` solves FHN traveling-wave equation and computes critical eigenfunction of the weighted eigenvalue problem. Select parameters and initial guess in lines 19-47 of this script and run this script in Matlab to reproduce the results in the paper cited above.

`FHNeqn_w_periodic.m` is an auxiliary script that returns the right-hand side of the FHN equation and the (weighted) Jacobian 

`fourdif.m` provides Fourier differentiation matrices

## Data Files

The initial data for the three different computations are contained in the Matlab data files `initial_pulse_{1,2,3}.mat`.
