These are the scripts used for the custom radar analyses underlying Culberg, et al "Double ridge formation over shallow 
water sills". 

Dependencies:
This code has been run using MATLAB R2019b Update 7 (9.7.0.1371314) running on 64-bit Windows 10 Pro Version 1909 and MATLAB
R2017b (9.3.0.713579) running on 64-bit CentOS linux release 7.7.1908 (Core). We expect this code should be compatible with 
any MATLAB release from R2017b onwards.

Installation:
Unzip to the desired directory. Installation time on a typical desktop computer is on the order of a few seconds. 

Directory Structure:
./ - top level scripts described below which peform the main analyses 
./DemoData - data used by the top level scripts for demonstration purposes. The scripts are configured to run directly on this 
             data.
./DerivedData - output data folder where results from scripts will be automatically saved
./RefData - reference data (such as climate model outputs) referenced by some of the scripts
./ReferenceFunctions - secondary functions used by some of the main level scripts to run their analyses
./SensitivityTests - scripts required to run the radar sensitivity tests described in Supplementary Methods - S5

Descriptions of Top Level Scripts:

Supplementary Methods - S2 - Subsurface Scattering Properties

RadarPhaseCoherence.m
	- This interactive scripts calculates the horizontal phase coherence of pulse compressed radar data at varying length scales and 
	  finds the mean coherence of different material types (ice, firn, and noise/clutter). To run interactively, place cursor in the 
	  portion of code to execute and hit CTRL-ENTER.
	- Time to Run: ~ 135s per radargram
	- Inputs and Instructions:
		- Download the pulse compressed radargram (PC_20160510_01_037.nc) from the Box link 
		  (https://stanford.box.com/s/2tip5unbdveddz7asdg0etmielgnbp5x) and save it in the ./DemoData directory before running the 
		  script.
	- Outputs:
		- phase coherence image
		- descriptive statistics for phase coherence of different material types in predefined regions of image
	- Data dependencies:
		- ./RefData/PC_20160510_01_037.mat
		- ./RefData/IceCoherenceBounds.mat
		- ./RefData/FirnCoherenceBounds.mat
		- ./RefData/NoiseCoherenceBounds.mat

		
Supplementary Methods - S3 - Absolute Calibration of Radar Data

The code used for the absolute calibration of the radar data has been previously published as part of:
R. Culberg, D. M. Schroeder, W. Chu, Extreme Melt Season Ice Layers Reduce Firn Permeability Across Greenland, 
Nature Communications, 2021. 10.1038/s41467-021-22656-5
The code is available at https://github.com/rtculberg/2012GreenlandMeltLayer and archived at https://zenodo.org/record/4552834.


Methods - Reflectivity Inversion, Supplementary Methods - S4 - Radar Forward Model 

CalculateAttenuationProfiles.main
	- This interactive scripts calculates the bounding attenuation profiles used to correct the radargrams based on firn core 
	  conductivity profiles and modeled temperature profiles.
	- Time to Run: ~ 2 s
	- Input:
		- All inputs are saved in the ./RefData directory and will be accessed automatically by the script.
	- Output:
		- d - depth axis from 0 to 50 meters (vector of doubles)
		- atten1 - total attenuation as a function of depth based on the MAR subsurface temperature profile (vector of doubles)
		- atten2 - total attenuation as a function of depth assuming the whole ice column is at the mean annual temperature (vector of doubles)
		- atten3 - total attenuation as a function of depth assuming that the whole ice column is at 0 C (vector of doubles)
	- Data dependencies:
		- ./RefData/B26_conductivity.mat
		- ./RefData/B27_conductivity.mat
		- ./RefData/B28_conductivity.mat
		- ./RefData/B29_conductivity.mat
		- ./RefData/MARTemperatureProfile.mat
		- ./RefData/MAR_20km_MeanVals_1980-2011.mat
			   
CalculateReflectivityDistributions.m
	- This interactive script calculates the mean and standard deviation of the reflectivity within predefined regions of a radargram.
	  The outputs of this script can be used as inputs to the MCMC reflectivity inversion.
	- Time to Run: ~69 s
	- Inputs and Instructions:
		- Download the SAR focused radargram (SAR_20160510_01_037.nc) from the Box link 
		  (https://stanford.box.com/s/2tip5unbdveddz7asdg0etmielgnbp5x) and save it in the ./DemoData directory before running the 
		  script.
	- Outputs:
		- Prints out the descriptive statistics for each predefined region of interest (mean, standard deviation, kurtosis, and skewness)
	- Data dependencies:
		- ./DemoData/SAR_20160510_01_037.nc
		- ./RefData/FirnBoundaries.mat
		- ./RefData/FirnBoundaries2.mat
		- ./RefData/IceBoundaries.mat
		- ./RefData/IceBoundaries2.mat
		- ./RefData/PlugIceBoundaries.mat

RunFacetMCMC.m
	- This function uses Markov Chain Monte Carlo methods with Metropolis-Hastings sampling to calculate the distributions subsurface
	  structure parameters consistent with observed radar reflectivity.
	- Time to Run: ~ 17 hrs for 1,000,000 iterations
	- Inputs:
		*** Note that all of these input values are available in ./DemoData/MCMCInputs.mat ***
		- data - mean radar reflectivity in dB (value to be inverted) (scalar double)
		- sigma - standard deviation of radar reflectivity distribution in dB (uncertainty in value to be inverted) (scalar double)
		- Niter - number of iterations for MCMC run (scalar double)
		- filename - name of output file to save results (string)
	- Outputs:
		- x_keep - 4 x Niter matrix of accepted model parameters
		- L_keep - 1 x Niter vector of the likelihood of the accepted model	parameters at that step
		- count - number of accepted proposals (count/Niter = acceptance rate)
		- reflectivity - 1 x Niter vector of modeled reflectivity from accepted parameters (can be compared to data distributions 
		  to sanity check the model runs)
	- Data dependencies:
		- The precalculated mean and standard deviations of the reflectivity for the firn, ice shell, and sub-ridge material are 
		  saved under ./DemoData/MCMCInputs.mat. 

MCMCResultsAnalysis.m
	- This interactive scripts takes the results of an MCMC inversion run and calculates the fractional ice content from the output 
	  parameters, then plots the histograms of those results.
	- Time to Run: ~2 s
	- Inputs:
		- Output of RunFacetMCMC.m - sample data is available in the ./DemoData directory.
	- Outputs:
		- Histograms and descriptive statistics for fractional ice content of each region of interest based on MCCM inversions.
	- Data Dependencies:
		- ./DemoData/MCMC_Firn_Final.mat
		- ./DemoData/MCMC_Ice_Final.mat
		- ./DemoData/MCMC_Ridge_Final.mat
		- ./DemoData/MCMC_Firn_SwapFirn.mat
		- ./DemoData/MCMC_Ice_SwapFirn.mat
		- ./DemoData/MCMC_Ridge_SwapFirn.mat


Supplementary Methods - S5 - Radar Inversion Sensitivity Tests

	- The ./SensitivityTests directory contains updated MCMC code for each of the sensitivity tests described in the supplementary 
	  methods. To run the sensitivity tests, the RunFacetMCMC.m function should be adjusted to call the appropriate version of 
	  mcmc_gauss_*.m as described below:
		- mcmc_gauss_fresnel.m - runs the MCMC inversion assuming that the scattering area is set by the diameter of the first Fresnel 
		  zone
		- mcmc_gaus_phaseerror.m - runs the MCMC inversion with a random phase screen that accounts for local surface roughness. When 
		  calling this version of the function, you will need to pass in addition variable "sigma_h" which is the local rms surface 
		  roughness in meters.
		- mcmc_gaus_swapfirn.m - runs the MCMC inversion assuming that the background material is firn with embedded patches of firn
		  and ice

