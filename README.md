# Updates

 * Version 1.1 (19June2017): Improved Eigenvalues/vectors calculations. Fixed other small bugs in `<leave_out_COMPLETE.m>`


# LeaveOutTwoWay
Leave Out estimates in two fixed effects models as described in Kline, Saggio and Sølvsten (2018).

Within this repository, there are three indepedent routines that users can test on their own
datasets:

Users interested in applying the homoskedastic correction in two-way
fixed effects models can use the function `<andrews_complete.m>`. 

Users interested in applying leave-out estimates in two-way
fixed effects models can use the function `<leave_out_COMPLETE.m>`.  

Users interested in computing leave-out estimates and inference 
of the variance of firm effects only can use the function 
`<leave_out_FD.m>` which works well with large datasets.

All is required for these functions to run appropriately is that the
original person-year file is appropriately sorted by worker-identifiers
and year.

Check the m-file `<main.m>` to run Leave Out estimates in a test sample. 

# Replication Archive
Link to replication archive of KSS: https://www.dropbox.com/sh/ozqq8qp9run9pzw/AAD1L9PZ-D4ueuIghq6eDh1Ea?dl=1 (last updated 06/20/2017)

