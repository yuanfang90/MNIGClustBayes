# MNIGClustBayes
A Bayesian approach for clustering skewed data using mixtures of multivariate normal-inverse Gaussian distributions.

A clustering algorithm based on mixture of MNIG distribution is proposed. Clustering and parameter estimation are done through a full Bayesian approach via Gibbs Sampling.

The 'program_2var.r' file provides the code to one simulation study where 2 groups of 2-variate MNIG observation are simulated and sent to the clustering algorithm. 'Simulation_2var_seed2019.RData' saves all results of one run of the code on the data simulated in 'program_2var.r'.

The 'program_4var.r' file provides the code to a simulation study where 3 groups of 4-variate MNIG observation are simulated and sent to the algorithm.
