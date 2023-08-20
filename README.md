# bim
A Born Iterative Method (BIM) implementation for use in Microwave Tomographic Imaging (MTI).

With this algorithm, the E-field results of COMSOL simulations are processed and the permittivity of each pixel is reconstructed to determine the presence of an anomaly.

The ultimate goal is for brain images to be reconstructed in terms of their dielectric properties to ascertain the presence and extent of brain trauma due to stroke or other factors.

## Setup
Cloning with Git: git clone https://github.com/hjbaum/bim.git

If using a ZIP, extract to your preferred folder.

## Use
(1) Ensure all required COMSOL simulation data has been exported to a csv and saved in a convenient location. 

Example: Simulation data for source 2 is located in \bim\COMSOL\Homogeneous\pt2_0_9GHz.csv

(2) Add file locations and their corresponding antenna IDs to filenames_by_source.csv located in \bim\Python\.

Example: The fields for source 2 can be written as 2	| 2	| Homogeneous	| pt2_0_9GHz_bleed_scatter.csv | pt2_0_9GHz_bleed_field.csv	| pt2_0_9GHz.csv | pt2_0_9GHz_field.csv	 

(3) Adjust the desired maximum # of iterations (MAX_ITER, line 208, bim.py), regularization factor (gamma, line 211, bim.py), operation frequency (freq, line 12, test.py) and other variables as needed. 

(4) Run test.py in the command prompt. 

(5) Wait for results to be plotted. 

## Main Resource
This algorithm is based on a 1989 paper by Chew and Wang.

The outline of the approach is as follows:

(1) Solve the linearized inverse problem for the first order distribution function by using the Born approximation. 

(2) Solve the scattering problem for the field in the object and at the observation points with the last reconstructed distribution function. 

(3) Substitute the field in the object obtained in step (2) into the integrand in the integral equation and solve the inverse problem to recover the next order distribution function. 

(4) Repeat step (2) and compare the field obtained by the reconstructed distribution function and the measured data, which in our case are the simulated fields for the exact distribution function at the observation points. If the difference is less than 5% of the scattered field, the iteration terminates, otherwise repeat the cycle until the solution is convergent.

Citation: Wang, Y. M., and W. C. Chew. “An iterative solution of the two-dimensional electromagnetic inverse scattering problem.” International Journal of Imaging Systems and Technology, vol. 1, no. 1, 1989, pp. 100–108, https://doi.org/10.1002/ima.1850010111. 

## TODO: 
### Rework for multiple sources (antennas). 
This has previously been attempted in branch wip_multiple_sources. This attempt was unsuccessful.

Main points:

(a) Call the BIM once from the test script using data collected from each antenna.

(b) K, scattered field, and epsilon (solution) data should be retained and refined with each iteration within the BIM.

(c) We should see a general improvement of the error and the solution as we perform calculations for each antenna's data, and throughout successive iterations.


### Continue debugging to determine sources of error, especially those based on a fundamental misunderstanding of the algorithm.
The current error is unreasonably and unexplicably high, likely due to the latter.

### Make the program modular and extensible
(a) Allow variables like the frequency, MAX_ITER, N, dx/dy, etc. to be set directly by the user via command prompt when running the program.

(b) Currently, test.py supports plotting the results of datasets from just 5 sources. The results for each source are overlaid in the final plot, and each set of results is assigned one of five colors defined in test.py line 122. Adding more colors to this list would extend the number of sources that can be plotted without causing an out-of-bounds exception. An alternative workaround would be to find a way of assigning colors dynamically and automatically with a python method. 

### Perform general cleanup
(a) Simplify redundant code. E.g. solve_baseline and solve_bleed in test.py are very similar in behavior and could be replaced by a single function.

(b) Reduce use of globals in bim.py where they are not needed.
