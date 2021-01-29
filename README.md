# StochasticParareal

This repository contains sample code for the pre-print paper by Pentland, Tamborrino, Sammadar & Appel - "Stochastic parareal: a novel application of probabilistic methods to time-parallelisation". 

This code is written in MATLAB and requires the use of the following toolboxes (noting there may be dependencies I have overlooked):
* Parallel Computing Toolbox.
* Statistics and Machine Learning.

## Authors

* Kamran Pentland - Mathematics Institute, University of Warwick
* Massimiliano Tamborrino - Dept. of Statistics, University of Warwick
* Debasmita Samaddar - Centre for Fusion Energy, Culham
* Lynton Appel - Centre for Fusion Energy, Culham


## Files
* Parareal: an implementation of the parareal algorithm.
* Stochastic Parareal: an implementation of the stochastic parareal algorithm.

Both files contain each algorithm respectively (parareal.m and stochasticparareal.m) as well as the RK.m function (an explicit Runge-Kutta solver) and TestProblems.m scripts (that implements each algorithm on various time-dependent ODE systems).



