# Brief Description

This repository computes Leave Out estimates of variance components in two fixed effects models as described in Kline, Saggio and Sølvsten (2018) - KSS henceforth - using Matlab. 

With the introduction of random projections techniques in version 2.0, it is possible to run leave out estimation of variance components on large datasets. To give an idea, for a dataset containing 5 million person year observations, 1.3 million person effects, 90K firm effects the code takes approximately 20 minutes to compute the relevant leave out matrices on a hpc system with 32 cores assigned.  

With the release of version 2.1, it is also possible to run inference on linear contrasts of regression coefficients that accounts for (i) many regressors asymptotics (ii) heteroskedasticity (iii) serial correlation within cluster. 

# Current Release: Version 2.1

It is now possible with "lincom_KSS" to make inference on linear combinations of regression coefficients suitable for settings with many regressors, heteroscedastic errors potentially serially correlated within "clusters". For theoretical justification and further details, see Proposition 1 and Remark 9 of KSS.

The typical example where one would apply "lincom_KSS" is when in a prior step, the user has estimated a two-way fixed effects
model and obtained the associated, say, worker and firm fixed effects. In a second step, the user is interested in regressing
the firm effects onto a bunch of observable characteristics. "lincom_KSS" provides the valid standard errors associated with this
regression.

See the function `codes/lincom_KSS` and `example_testing` for a set of different scenarios where one would be interested in using the function "lincom_KSS".
 
# History of Updates

 * Version 2.0.1: Fixed Bug where code could not run because it could not locate the function `leave_out_estimation_two_way'.

 * Version 2.0: Added option to estimate the stastistical leverages, Pii, and what we define as Bii in KSS using Random Projections methods that build on the Johnson Lindestrauss Lemma - See Appendix B of KSS.
 
     * This especially helpful in massive datasets where exact computation of (Bii,Pii), even after the improvements introduced from version 1.5, is close to be prohibitive in terms of computation time.
 
     * A benchmark: in a data set with approximately 5 million observations, 1.3 million worker effects, and 90 thousand firm effects, it takes *21 minutes* to compute (Bii,Pii) via random projections. By contrast, it takes more than *30 hours* using the exact method." 
 
     * In terms of the code, We added the following inputs
                
                *   "type_of_algorithm": This takes two values: "exact or "JLL".

                           "type_of_algorithm = exact": performs exact computation of (Bii,Pii)
                           as described in version 1.5.          
                    
                           "type_of_algorithm = JLL": performs random projection methods to
                           approximate (Bii,Pii) as detailed in Appendix B of KSS.
                
               *    "epsilon": this governs the tradeoff b/w speed and accuracy 
                    when estimating (Bii,Pii). Smaller values of epsilon implies 
                    more accuracy but slower performance. See the paper by 
                    Spielman-Srivastava (2011) for further details.  

     * For a more complete picture of the accuracy vs. speed tradeoff involving the random projection approach, here are some diagnostics  for a dataset with approx 50K workers, 15K Firms

                *  When "type_of_algorithm = exact", the code takes 1500
                   seconds to compute (Bii,Pii) for variance of firm
                   effects, variance of person effects and covariance of
                   person, firm effects. 
 
                *  When "type_of_algorithm = JLL" and `epsilon`=0.005, 
                   the code takes 82 seconds to compute (Bii,Pii) 
                   for variance of firm effects, variance of person effects 
                   and covariance of person, firm effects. 
                   
                *  The correlation b/w the P_ii found with the exact method and 
                   the Pii found with the random projection method is equal
                   to 0.9987. The maxium difference found between the two is
                   less than <0.05.
 
                *  The variance components estimates from JLL differ from those 
                   obtained with the exact algorithm by a factor less than 0.1 percent.
                   Similarly, for the estimated standard errors (subsample_llr_fit=0).

 * Version 1.55 (22August2018): Added more options to run non-parametric fit.
 
 * Version 1.52 (21August2018): Added more outputs to the function `leave_out_COMPLETE.m`  to simplify possible post-estimation commands.

 * Version 1.51 (15August2018): Better management of large sparse matrices when invoking parfor to compute (Bii,Pii) using the option      `parallel.pool.Constant`.

 * Version 1.5 (14August2018): Improved "exact" computation of (Bii,Pii). In particular, we have introduced the following changes:
                
                * Added CMG routine to speed computation of linear system
                  involving the Laplacian matrix as design matrix. 
                
                * CMG package - available here: http://www.cs.cmu.edu/~jkoutis/cmg.html - 
                  has been already included in the repository.
 
                * Read movers-stayers structure to fasten computation of (Bii,Pii).
                
  * In terms of speed, for the test dataset used in `example.m`
  
                * With version 1.32 the code takes 260 seconds to compute (Bii,Pii).
                
                * With version 1.5 the code takes 23 seconds to compute (Bii,Pii).
                
 * Version 1.32 (01JAug2018): Introduced example.m with better management of folders where results are saved. `leave_out_COMPLETE.m` now                               also exports a .csv file containing the main variables in the leave out connected set
 * Version 1.31 (31Jul2018): Added the option "DO_SE" in main.m.
 * Version 1.3 (25Jul2018): Dropped stayers with a single person-year observations (for whom Pii=1 when estimating model in levels).
 * Version 1.2 (01Jul2018): Added more options to speed-up computation of the standard errors.
 * Version 1.13 (22June2018):Better read of Nargout options. Added option if user wants computation of standard error.
 * Version 1.12 (20June2018): Fixed minor bugs when running leave out estimation with controls. 
 * Version 1.1 (19June2018): Improved Eigenvalues/vectors calculations. Fixed other small bugs in `leave_out_COMPLETE.m`
 * Version 1.0 (16June2018): Original upload. 

# Routines included within the Repository.
Within this repository, there are four independent routines that users can test on their own
datasets:

* Users interested in applying the homoskedastic correction in two-way
fixed effects models can use the function `andrews_complete.m`. 

* Users interested in applying leave-out estimates in two-way
fixed effects models can use the function `leave_out_COMPLETE.m`.  

* Users interested in computing leave-out estimates and inference 
of the variance of firm effects only can use the function 
`leave_out_FD.m`.

* Users interested in computing valid inference on linear combination 
of regression coefficients that is robust to many regressors asymptotics, 
heteroskedasticity and serial correlation within clusters can use
`lincom_KSS.m`. lincom_KSS can be applied to ANY linear regression model (i.e.
not only two-way models) but it optimized to run as essentially a post-estimation
command following the command `leave_out_COMPLETE.m`.

# Things to have in mind 

* `example.m` provides a simple example where the user calls a test dataset in .csv and then calls the function `leave_out_COMPLETE.m` to compute leave out estimates. See the documentation provided within the function of `leave_out_COMPLETE.m` for a description of the inputs and the outputs associated with this function.

* `example_testing.m` provides five examples that show how to conduct inference on linear combination of linear regression parameters using the function `lincom_KSS.m.

* `leave_out_COMPLETE.m` outputs estimates of the VCM of the person, firm effects with associated standard errors (assuming q=0 in KSS).  Standard errors for the case q=1 can also be obtained. In future releases, we plan also to allow users to specificy the variance component associated with a specific (sub)set of controls. 

* Adding Controls to the model: This is easily handled by `leave_out_COMPLETE.m`. However, to speed up computations, we suggest to partial out the effects of these controls in a first step so that the associated design matrix has a Laplacian structure. See the documentation provided inside the functions `leave_out_COMPLETE.m` and `leave_out_FD.m` and Appendix B of KSS.

* It is highly suggested that the user runs any of the function mentioned above on an hpc system in order to fully exploit the parallelization features of Matlab.

* Code runs on any type of person-year file (balanced vs. unbalanced; T=2 vs. T>2). The only important thing is that the variables in memory (outcome variable, worker identifiers, firm identifiers, eventual controls) need to be sorted by worker id and year.

# Replication Archive
Link to replication archive of KSS: https://www.dropbox.com/sh/ozqq8qp9run9pzw/AAD1L9PZ-D4ueuIghq6eDh1Ea?dl=1 (last updated 06/20/2017)

