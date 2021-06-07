# TRINITY
An empirical model that statistically connects dark matter halos, galaxies and supermassive black holes (SMBHs) from z=0-10.

## Prerequisite
This project requires OpenMP and GSL. The users should specify the path to GSL by setting up environmental variables $GSL_INCLUDE and $GSL_LIB:

$GSL_INCLUDE=/your/path/to/gsl/include

$GSL_LIB=/your/path/to/gsl/lib

## Installation
The user can use the Makefile to compile the whole project. (TBC)

## Roadmap
After successful compilation, multiple binary executables will be produced:

1. fitter is used for finding the best fitting parameter using a iterative gradient descent algorithm;
2. all_smf_mcmc is used for performing Markov Chain Monte Carlo;
3. Other files, e.g., gen_smf and gen_edd_r_color, are used for predicting galaxy and SMBH properties based on input model parameters.

## Renerating Predictions Based on Model Parameters
TRINITY is able to are predict many observational data (e.g., galaxy stellar mass functions and quasar luminosity functions) and underlying galaxy and SMBH properties (e.g., SMBH Eddington average Eddington ratios). These predictions are made by different code files like gen_smf.c, gen_qlf.c, gen_edd_r_color.c, etc.. There are broadly two types of ''prediction'' codes: the first type generate observable data given input redshift or redshift invertals like gen_smf.c and gen_qlf.c; the second type generates galaxy or SMBH properties as a function of host halo mass and redshift, like gen_edd_r_color.c.

### Predictions of observables like stellar mass functions and luminosity functions (Figs. 3-9 of Zhang et al. 2021)
Below is how to generate quasar luminosity functions with gen_qlf:

./gen_qlf z mass_cache mcmc_output > qlf.dat

Here, z is the input redshift at which the quasar luminosity function will be calculated, and mass_cache is the file containing the cached halo mass functions (we used mf_bolshoi_planck.dat in TRINITY/aux/ for Zhang et al. 2021). mcmc_output is the string containing the model parameters, which can be found in TRINITY/aux/best_fit_final_eff_z.dat. qlf.dat is the output file containing the predicted quasar luminosity function.


### Predictions of observables like stellar mass functions and luminosity functions (Figs. 15-19, 23 of Zhang et al. 2021)
Below is how to predict galaxy and SMBH properties as functions of halo mass and redshift with codes with gen_edd_r_color:

./gen_edd_r_color mass_cache (mcmc output) > properties.dat.

The required inputs are similar to the usage of gen_qlf except for the redshift.
