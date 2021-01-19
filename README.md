# LeaveOutTwoWay

This Matlab package implements the leave out correction of 
[Kline, Saggio and Soelvsten (2020)](https://eml.berkeley.edu/~pkline/papers/KSS2020.pdf) for estimating variance components in two-way fixed effects models. 
 
See this [vignette](doc/VIGNETTE.pdf) for a description of the package.

## Julia Version

The Julia version of the package can be 
found [here](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl). On the same
page, one can find an executable that permits estimation of the leave-out correction even if the user does not have MATLAB or JULIA. 

## Summary of Changes

These are the most significant changes introduced with this new version
* New documentation that describes in detail the functioning of `leave_out_KSS`.
* By default, the code runs a leave-out correction by leaving a match out as opposed to leaving an observation out. See [vignette](doc/VIGNETTE.pdf) for details.
* New, optimized, random projection algorithm for calculation of the statistical leverages which scale extremely well to large datasets.
* Code no longer requires MATLAB BGL but runs automatically on MATLAB built-in network functions.
* Requires MATLAB R2015b or higher. 

## Replication Package
The replication package of the Econometrica article can be found [here](https://www.dropbox.com/s/iaj3aap3ibfhup8/Replication%20ECTA.zip?dl=1). 



