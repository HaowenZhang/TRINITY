# TRINITY
TRINITY is an empirical model that statistically connects dark matter halos, galaxies and supermassive black holes (SMBHs) from z=0-10. Constrained by multiple galaxy (0<z<10) and SMBH datasets (0<z<6.5), TRINITY finds the posterior probability distributions of the halo--galaxy--SMBH connection and SMBH properties, all of which are allowed to evolve with redshift.

Most code: Copyright (C)2018-2021 Haowen Zhang

License: GNU GPLv3

Science/Documentation Paper: https://arxiv.org/abs/2105.10474

## Prerequisites
This project requires OpenMP and GSL. The users should specify the path to GSL by setting up environmental variables $GSL_INCLUDE and $GSL_LIB:

$GSL_INCLUDE=/your/path/to/gsl/include

$GSL_LIB=/your/path/to/gsl/lib

## Installation
TRINITY comes with a Makefile, so user can use the "make" command to compile the whole project. Alternatively, users can use target names specified in the Makefile to specifically produce the corresponding executables. For example, if one wants to compile the code to predict quasar luminosity functions (QLFs) only, there is a target in the Makefile as follows:


qlf:
	
$(CC) $(BASE_FILES) gen_qlf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf
  
So they can use the command "make qlf" to compile gen_qlf.c and produce gen_qlf, without having to compile the entire project.

## Roadmap
After successful compilation, multiple binary executables will be produced:

1. TRINITY/src/fitter is used for finding the best fitting parameter using a iterative gradient descent algorithm;
2. TRINITY/src/all_smf_mcmc is used for performing Markov Chain Monte Carlo (MCMC);
3. Other files, e.g., gen_smf and gen_edd_r_color, are used for predicting galaxy and SMBH properties based on input model parameters. For the purpose of each executables, please see the comments in their corresponding c codes.

## Renerating Predictions Based on Model Parameters
TRINITY is able to are predict many observational data (e.g., galaxy stellar mass functions and quasar luminosity functions) and underlying galaxy and SMBH properties (e.g., SMBH Eddington average Eddington ratios). These predictions are made by different code files like gen_smf.c, gen_qlf.c, gen_edd_r_color.c, etc.. There are broadly two types of ''prediction'' codes: the first type generate observable data given input redshift or redshift invertals like gen_smf.c and gen_qlf.c; the second type generates galaxy or SMBH properties as a function of host halo mass and redshift, like gen_edd_r_color.c.

### Predictions of Observables Such As Galaxy Stellar Mass Functions and Quasar Luminosity Functions (Figs. 3-9 of Zhang et al. 2021)
Below is how to generate quasar luminosity functions with gen_qlf:

./gen_qlf z mass_cache param_file > qlf.dat

Here, z is the input redshift at which the quasar luminosity function will be calculated, and mass_cache is the file containing the cached halo mass functions (we used TRINITY/aux/mf_bolshoi_planck.dat in Zhang et al. 2021). param_file is the file containing the model parameters, e.g., TRINITY/aux/best_fit_final_eff_z.dat. qlf.dat is the output file containing the predicted quasar luminosity function.

### Predictions of Galaxy and SMBH Properties as Functions of Halo Mass and Redshift (Figs. 15-19, 23 of Zhang et al. 2021)
Below is how to predict galaxy and SMBH properties as functions of halo mass and redshift with codes with gen_edd_r_color:

./gen_edd_r_color mass_cache param_file > properties.dat.

The required inputs are similar to the usage of gen_qlf except for the redshift.

## Running Your Own Model Fitting or Markov Chain Monte Carlo (MCMC) Processes
./fitter carries out model fitting with gradient descent method. The usage is shown as follows (TRINITY/examples/fit.pbs):

OMP_NUM_THREADS=8 ../src/fitter 1 1 1 1 1 0 0 ../aux/mf_bolshoi_planck.dat  ../obs/Kelly_BHMF_TypeI/\*.bhmf_typei ../obs/\*.\* ../obs/Ueda_QLF_CTK/\*.qlf ../obs/Aird_qpdf_no_high/\*no_highest_0.3dex.qpdf ../obs/qf/\*.qf ../obs/uvlf/new/\*.uvlf < ./cfit_final_eff_z.param > ../mcmc_runs/fit_final_eff_z.dat 2> ../mcmc_runs/err_final_eff_z.dat

Here, several numerical flags are set (1 1 1 1 1 0 0). For the meaning of each flag, see init_mcmc_from_args() in all_smf.c. ./fitter takes observational data (e.g., those included in TRINITY/obs/) as constraints, and read initial parameters from stdin. The final output parameters are directed to stdout and the log messages to stderr.

./all_smf_mcmc runs MCMC process using adaptive Metropolis-Hastings method. The usage is shown as follows (TRINITY/examples/mcmc.pbs): 

OMP_NUM_THREADS=16 ../src/all_smf_mcmc 1 1 1 1 1 0 0 ../aux/mf_bolshoi_planck.dat ../obs/Kelly_BHMF_TypeI/\*.bhmf_typei ../obs/\*.\* ../obs/Ueda_QLF_CTK/\*.qlf ../obs/Aird_qpdf_no_high/\*no_highest_0.3dex.qpdf ../obs/qf/\*.qf ../obs/uvlf/new/\*.uvlf < ../mcmc_runs/fit_final_eff_z.dat > ../mcmc_runs/all_final_eff_z.dat 2> ../mcmc_runs/burn_in_final_eff_z.dat

The format is very similar to the model fitting case.
