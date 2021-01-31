# StochasticParareal

This repository contains sample code for the pre-print paper by Pentland, Tamborrino, Sammadar & Appel - "Stochastic parareal: a novel application of probabilistic methods to time-parallelisation". 

This code is written in MATLAB and requires the use of the following toolboxes (noting there may be dependencies I have overlooked):
* Parallel Computing Toolbox.
* Statistics and Machine Learning.

## Authors

* Kamran Pentland - Mathematics Institute, University of Warwick
* Massimiliano Tamborrino - Dept. of Statistics, University of Warwick
* Debasmita Samaddar - Cuham Centre for Fusion Energy, Abingdon, Oxfordshire
* Lynton Appel - Cuham Centre for Fusion Energy, Abingdon, Oxfordshire


## Files
* Parareal: an implementation of the parareal algorithm.
* Stochastic Parareal: an implementation of the stochastic parareal algorithm.

Both files contain each algorithm respectively (parareal.m and stochasticparareal.m) as well as the RK.m function (an explicit Runge-Kutta solver that both main algorithms use) and TestProblems.m scripts (scripts that give example use cases for each algorithm on various time-dependent ODE systems).

To run the algorithms, open either of the TestProblems.m scripts and run one of the sections of code containing a particular test problem. Each section plots the solutions to the ODEs using both the Rk.m solver and either parareal.m or stochasticparareal.m. 



