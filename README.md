# Updates

 * Version 1.1 (19June2017): Improved Eigenvalues/vectors calculations. Fixed other small bugs in `<leave_out_COMPLETE.m>`
 * Version 1.0 (16June2017): Original upload. 


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
and year (xtset id year in Stat).

The m-file `main.m` runs Leave Out estimates in an test dataset provided within this repository. 

# Things to have in mind 

* `leave_out_COMPLETE.m` outputs estimates of the VCM of the person, firm effects with associated standard errors (assuming q=0 in KSS).  Standard errors for the case q=1 can also be obtained. In future releases, we plan also to allow users to specificy the variance component associated with a specific (sub)set of controls. 

* `leave_out_COMPLETE.m` is currently not optimized for large datasets. We'll make an effort to constantly update the code in order to deal with very large datasets. Version 1.1 is a first step towards that direction. 

* Adding Controls to the model: This is easily handled by `leave_out_COMPLETE.m`. However, to speed up computations, we suggest to partial out the effects of these controls in a first step so that the associated design matrix has a Laplacian structure. See the documentation provided inside the functions `leave_out_COMPLETE.m` and `leave_out_FD.m` and Appendix B of KSS.

* Code runs on any type of person-year file (balancend-unbalanced), the only important thing is that the variables in memory (outcome variable, worker identifiers, firm identifiers, eventual controls) need to be sorted by worker id and year.

# Replication Archive
Link to replication archive of KSS: https://www.dropbox.com/sh/ozqq8qp9run9pzw/AAD1L9PZ-D4ueuIghq6eDh1Ea?dl=1 (last updated 06/20/2017)

