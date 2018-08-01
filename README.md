# Next Releases

 * Next release will allow for fast computation of covariance of person, firm effects and variance of person effects in large datasets. 

# Updates

 * Version 1.31 (31Jul2018): Added the option "DO_SE" in main.m.
 * Version 1.3 (25Jul2018): Dropped stayers with a single person-year observations (for whom Pii=1 when estimating model in levels).
 * Version 1.2 (01Jul2018): Added more options to speed-up computation of the standard errors.
 * Version 1.13 (22June2018):Better read of Nargout options. Added option if user wants computation of standard error.
 * Version 1.12 (20June2018): Fixed minor bugs when running leave out estimation with controls. 
 * Version 1.1 (19June2018): Improved Eigenvalues/vectors calculations. Fixed other small bugs in `leave_out_COMPLETE.m`
 * Version 1.0 (16June2018): Original upload. 


# Description 
Leave Out estimates in two fixed effects models as described in Kline, Saggio and Sølvsten (2018) - KSS henceforth.

Within this repository, there are three independent  routines that users can test on their own
datasets:

* Users interested in applying the homoskedastic correction in two-way
fixed effects models can use the function `andrews_complete.m`. 

* Users interested in applying leave-out estimates in two-way
fixed effects models can use the function `leave_out_COMPLETE.m`.  

* Users interested in computing leave-out estimates and inference 
of the variance of firm effects only can use the function 
`leave_out_FD.m` which works well with large datasets.

All that is required for these functions to run appropriately is that the
original person-year file is sorted by worker-identifiers
and year (xtset id year in Stata).

The m-file `example.m` runs Leave Out estimates in an test dataset provided within this repository. 

# Things to have in mind 

* `example.m` provides a simple example where the user calls a test dataset in .csv and then calls the function `leave_out_COMPLETE.m` to compute leave out estimates. See the documentation provided within the function of `leave_out_COMPLETE.m` for a description of the inputs and the outputs associated with this function.

* `leave_out_COMPLETE.m` outputs estimates of the VCM of the person, firm effects with associated standard errors (assuming q=0 in KSS).  Standard errors for the case q=1 can also be obtained. In future releases, we plan also to allow users to specificy the variance component associated with a specific (sub)set of controls. 

* Adding Controls to the model: This is easily handled by `leave_out_COMPLETE.m`. However, to speed up computations, we suggest to partial out the effects of these controls in a first step so that the associated design matrix has a Laplacian structure. See the documentation provided inside the functions `leave_out_COMPLETE.m` and `leave_out_FD.m` and Appendix B of KSS.

* It is highly suggested that the user runs any of the function mentioned above on an hpc system in order to fully exploit the parallelization features of Matlab.

* Code runs on any type of person-year file (balancend-unbalanced), the only important thing is that the variables in memory (outcome variable, worker identifiers, firm identifiers, eventual controls) need to be sorted by worker id and year.

# Replication Archive
Link to replication archive of KSS: https://www.dropbox.com/sh/ozqq8qp9run9pzw/AAD1L9PZ-D4ueuIghq6eDh1Ea?dl=1 (last updated 06/20/2017)

