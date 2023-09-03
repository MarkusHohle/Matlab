# Matlab
Pattern formation, diffusion reaction, Bayesian 

Set of Matlab codes from various lectures

# Bayesian Statistics

- bayesian_bino.m (Bayesian Parameter Estimation of a Binomial Distribution)
- bayesian_poiss.m (Bayesian Parameter Estimation of a Poisson Distribution)

- trapezoidal_rule_nd_integral.m (function for N-dimensional numerical integration. Author: Mohammed S. Al-Rawi, IEETA, University of   Aveiro Portugal)

- model_selection_example.m (performs Model Selection on 'tumor animal Fig2.txt', data taken from DOI:10.1093/carcin/23.3.511, calls    trapezoidal_rule_nd_integral.m)

- mean_var_analysis_example.m (performs Model Selection on 'LinaData.mat': WT vs treated mice, returns odds ratio and compares result   to standat t-test. Data with permission from Dr. Lina Wendeler)

- VarBayesExample.m (estimating mean and variance and calculates corresponding probability density distributions using Variational      Bayes, calls same data set as above)

- jumpGraph.m & jumpGraphRandom.m (illustrates temporal changes of states on a graph network with fixed rates and randomly generated    rates which get multiplied with the adjacency matrix; calls axesArea.m, prh.m and wgPlot.m for visualisation. Last three functions    have been created by Michael Wu waftingpetal@yahoo.com)


# Biophysics - ODEs for non-linear systems

- poiss_stepper.m (simulates poissonian stepper with a for loop)
  
- poiss_stepper_fast.m (simulates poissonian stepper with matrix multiplication)

- poiss_stepper_manyTimes.m (simulates N poissonian stepper in order to illustrate the law of large numbers)

- random_machine.m (simulates the inevitable increase of entropy as a consequence being the most likely macrostate with N dice          getting rolled M times)

- random_machine_Boltzmann.m (same as above, but with the constrain energy = const  --> Boltzman distribution emerges, see Lagrangian   Multipliers)

- glycolysis.m (ODEs for simplyfied glycolysis reaction), solve_glycolysis.m (solves corresponding ODS and plots result, i. e. fixed    points etc, calls plotpowerspec.m )

- glycolysis_2.m, solve_glycolysis_2.m (same as above, but reaction rates a and b can be provided as input arguments)

- plotpowerspec.m (performs FFT and plots frequency power spectrum)

- phage.m, solve_phage.m (simulates immune reaction/bacteria growth as response of phage therapy, creates plots, see the work from      DOI: 10.1016/j.jtbi.2017.06.037)

- gillespieI.m - gillespieIV.m (simulates different systems like predator - prey, glycolysis etc as stochastic process using            Gillespie algorithm)

- example_gene_expression_gillespie.m (simulates gene expression, as stochastic process, author: nezar@mit.edu, central dogma,          master equation: see Physics of Life Reviews 2 (2005) 157–175, equation 3, calls firstReactionMethod.m)

- thattai_ode.m (solves deterministic ODEs of mRNA - protein interaction, Author: Erwin Frey, frey@lmu.de)



# Biophysics - Stochastic Processes

- random_walk.m (simple illustration of a random walk)

- biased_random_walk.m (simulates trajectories of N E-coli following a food concentration gradient)

- diffusion_reaction.m, diffusion_reactionII.m, diffusion_reactionIII.m (different examples for pattern formation like fur, scales      etc, see Koch & Meinhardt 1994)

- twoDdiffusion.m (generic example)

- twoDdiffusion_smoluchouski.m (diffusion with drift i. e. Smoluchouski equation)



