# TRINITY
TRINITY is an empirical model that statistically connects dark matter halos, galaxies and supermassive black holes (SMBHs) from z=0-10. Constrained by multiple galaxy (0<z<10) and SMBH datasets (0<z<6.5), TRINITY finds the posterior probability distributions of the halo--galaxy--SMBH connection and SMBH properties, all of which are allowed to evolve with redshift.

Most code: Copyright (C)2018-2021 Haowen Zhang

License: GNU GPLv3

Science/Documentation Paper: https://arxiv.org/abs/2105.10474

## 1. Prerequisites
This project requires OpenMP and GSL. The users should specify the path to GSL by setting up environmental variables $GSL_INCLUDE and $GSL_LIB:

$GSL_INCLUDE=/your/path/to/gsl/include

$GSL_LIB=/your/path/to/gsl/lib

Alternatively, the users can also modify the $GSL_FLAG variable in Makefile to include these dependencies.

## 2. Installation
TRINITY comes with a Makefile (TRINITY/src/Makefile), so user can use the "make" command to compile the whole project. Alternatively, users can use target names specified in the Makefile to specifically produce the corresponding executables. For example, if one wants to compile the code to generate plots in Paper I, there is a target (i.e., named groups of codes to be compiled) called paper1_plots. To compile these codes only, users can simply execute the following command:
	
	make paper1_plots
  
. For the full list and content of all the compiling targets, please see TRINITY/src/Makefile.

## 3. Roadmap
1. TRINITY/src contains the source codes of TRINITY:
2. TRINITY/src/base contains the basic utility codes that are used in all fitting, MCMC, and prediction codes;
3. TRINITY/src/fitting_mcmc contains the codes to run gradient descent optimization of model parameters (the code to run MCMC is all_smf.c, which is in TRINITY/src/base);
4. TRINITY/src/paper1 contains the codes to generate the plots in Paper I;
5. TRINITY/src/paper2 contains the codes to generate the plots in Paper II;
6. The codes to generate plots for future papers will be added to the repository upon the acceptance of those papers.


After successful compilation, multiple binary executables will be produced:

1. TRINITY/bin/fitter is used for finding the best fitting parameter using a iterative gradient descent algorithm;
2. TRINITY/bin/all_smf_mcmc is used for performing Markov Chain Monte Carlo (MCMC);
3. Other files, e.g., gen_smf and gen_edd_r_color, are used for predicting galaxy and SMBH properties based on input model parameters. For the purpose of each executables, please see the comments in their corresponding c codes.

## Renerating Predictions Based on Model Parameters
TRINITY is able to are predict many observational data (e.g., galaxy stellar mass functions and quasar luminosity functions) and underlying galaxy and SMBH properties (e.g., SMBH Eddington average Eddington ratios). These predictions are made by different code files like gen_smf.c, gen_qlf.c, gen_edd_r_color.c, etc.. There are broadly two types of ''prediction'' codes: the first type generate observable data given input redshift or redshift invertals like gen_smf.c and gen_qlf.c; the second type generates galaxy or SMBH properties as a function of host halo mass and redshift, like gen_edd_r_color.c.

### Predictions of Observables Such As Galaxy Stellar Mass Functions and Quasar Luminosity Functions (Figs. 3-9 of Paper I)
Below is how to generate quasar luminosity functions with gen_qlf:

./gen_qlf z mass_cache > qlf.dat

Here, z is the input redshift at which the quasar luminosity function will be calculated, and mass_cache is the file containing the cached halo mass functions (we used TRINITY/aux/mf_bolshoi_planck.dat in Zhang et al. 2021). qlf.dat is the output file containing the predicted quasar luminosity function. Unless otherwise noted, the default model parameters used in such calculations are the best-fitting values of the fiducial TRINITY model (see Paper I and the erratum), which can be found in TRINITY/src/base/params_default.h.

### Predictions of Galaxy and SMBH Properties as Functions of Halo Mass and Redshift (Figs. 15-19, 23 of Paper I)
Below is how to predict galaxy and SMBH properties as functions of halo mass and redshift with codes with gen_edd_r_color:

./gen_edd_r_color mass_cache > properties.dat.

The required inputs are similar to the usage of gen_qlf except for the redshift.

## Running Your Own Model Fitting or Markov Chain Monte Carlo (MCMC) Processes
./fitter carries out model fitting with gradient descent method. The usage is shown as follows (TRINITY/examples/fit.pbs):

OMP_NUM_THREADS=8 ../src/fitter 1 1 1 1 1 0 0 ../aux/mf_bolshoi_planck.dat  ../obs/Kelly_BHMF_TypeI/\*.bhmf_typei ../obs/\*.\* ../obs/Ueda_QLF_CTK/\*.qlf ../obs/Aird_qpdf_no_high/\*no_highest_0.3dex.qpdf ../obs/qf/\*.qf ../obs/uvlf/new/\*.uvlf < ./cfit_final_eff_z.param > ../mcmc_runs/fit_final_eff_z.dat 2> ../mcmc_runs/err_final_eff_z.dat

Here, several numerical flags are set (1 1 1 1 1 0 0). For the meaning of each flag, see init_mcmc_from_args() in all_smf.c. ./fitter takes observational data (e.g., those included in TRINITY/obs/) as constraints, and read initial parameters from stdin. The final output parameters are directed to stdout and the log messages to stderr.

./all_smf_mcmc runs MCMC process using adaptive Metropolis-Hastings method. The usage is shown as follows (TRINITY/examples/mcmc.pbs): 

OMP_NUM_THREADS=16 ../src/all_smf_mcmc 1 1 1 1 1 0 0 ../aux/mf_bolshoi_planck.dat ../obs/Kelly_BHMF_TypeI/\*.bhmf_typei ../obs/\*.\* ../obs/Ueda_QLF_CTK/\*.qlf ../obs/Aird_qpdf_no_high/\*no_highest_0.3dex.qpdf ../obs/qf/\*.qf ../obs/uvlf/new/\*.uvlf < ../mcmc_runs/fit_final_eff_z.dat > ../mcmc_runs/all_final_eff_z.dat 2> ../mcmc_runs/burn_in_final_eff_z.dat

The format is very similar to the model fitting case.

## TRINITY Paper II: Predicting the AGN luminosity-dependent bias in the SMBH mass--galaxy mass relation

### 1. Compilation
In ./src/, users can use the command "make paper2_plots" to compile the codes to reproduce Figs. 1-2 of TRINITY Paper II (https://doi.org/10.1093/mnrasl/slad060). The compiled binary files are stored in ./bin/ . 

### 2. Luminosity-dependent SMBH mass--galaxy mass relation (BHSM) as a function of redshift and lower limit in AGN bolometric luminosity
The file ./bin/mbh_perc_mstar_lbol predicts the 16th, 50th, and 84th percentiles of SMBH mass as a function of galaxy stellar mass and redshift for AGNs/quasars above certain bolometric luminosities. The usage is as follows: 

	./bin/mbh_perc_mstar_lbol z log(lbol_min) sigma_mbh mf_cache >./output.dat
, where z is the redshift, log(lbol_min[erg/s]) is the log10 of the lower limit in AGN bolometric luminosity, and sigma_mbh is the random scatter in the log10 of observed SMBH mass around the intrinsic value. Typically, 0.5 dex is adopted for virial estimates. mf_cache is a file containing cached halo mass functions (HMFs) that are pre-calculated based on dark matter N-body simulations. We also included the cached HMF file we used in this repository: ./aux/mf_bolshoi_planck.dat.

Below is an example to generate the BHSM relation (including 16th, median, and 84th percentile values) for z=6 quasars brighter than $\log L_\mathrm{bol}>=46.5$, with a 0.5 dex random scatter in observed SMBH mass:

	./bin/mbh_perc_mstar_lbol 6.0 46.5 0.5 ./aux/mf_bolshoi_planck.dat >./bhsm_quasars.dat

 
### 3. Deviation from the quasar SMBH mass--galaxy mass relation
The file ./bin/mbh_perc_mstar_lbol_individual calculates the deviation in SMBH mass from the expected quasar BHSM relation, for individual quasars. The usage is as follows: 

	./bin/mbh_perc_mstar_lbol_individual mf_cache quasar_catalog >./output.dat
 , where quasar_catalog is a (text file, not binary) catalog of quasars with the following format:

 Column 1: redshift
 
 Column 2: log10 of SMBH mass [Msun]
 
 Column 3: log10 of galaxy mass [Msun]
 
 Column 4: log10 of AGN bolometric luminosity [erg/s]

 Column 5: uncertainty in SMBH mass [dex]

The first four(4) columns of output.dat will be the same as quasar_catalog.dat, but the fifth column will contain the deviation in SMBH mass from the expected quasar BHSM relation AT each quasar/AGN's bolometric luminosity. NOTE the difference from ./mbh_perc_mstar_lbol, which calculates BHSM relations ABOVE a certain luminosity. For the full details between these two calculations, we refer readers to 
	
 	./src/paper2/mbh_perc_mstar_lbol(_individual).c.

