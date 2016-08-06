For the algorithm description check "Approximate Bayesian Computation by Subset Simulation" by Manuel Chiachio, James L. Beck, Juan Chiachio, and Guillermo Rus, SIAM J. Sci. Comput., 36(3), A1339-A1358.

To work with the ABC-SubSim one needs to provide a fitfun.cpp file describing the simulation and edit the 'problem.h' file. 
Important things defined in these files are:

fitfun.cpp:
- experimental data used for calibration
- a 'get_discrepancy' function which compares the experimental data and the simulation outcomes

problem.h:
- parameter space specifications (dimensionality of the problem, bounds for the prior distribution, prior distribution specifications),
- 'run_problem' function which calls the simulation.

The code uses the GSL library and the TORC task-stealing library.

Once the code is compiled it can be run with './<execname>' to print the list of command line arguments one needs to provide.
The arguments are: <population size> <posterior filename> <final tolerance> <max #levels> <acceptance rate> <weights std> <tolerance change> <random seed> [<restart filename>].

- 'restart filename' argument is optional and provides the opportunity to restart the code using one of the output files produced by itself
  (the file must contain N = population size lines with a sample and a corresponding discrepancy  in each line).
- 'random seed' is meant for reproducibility.
- 'population size' is self-explanatory.

Other arguments are various stopping criteria. Check the generated output to see their values.

- 'final tolerance': stop when this tolerance reached
- 'max #levels': stop after that many generations
- 'acceptance rate': stop when MCMC acceptance rate drops below this threshold
- 'weights std': stop when weights standard deviation drops below this threshold
- 'tolerance change': stop when the relative change of tolerance drops below this threshold

The recommended stopping criterium is 'acceptance rate' (good value is ~0.05) which means that all the rest (except for 'max #levels') can be set to 0.

For details or help you can ask: lina.kulakova AT mavt.ethz.ch
