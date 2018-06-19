# LeaveOutTwoWay
Leave Out estimates in two fixed effects models as described in Kline, Saggio and Sølvsten (2018)

Here is a description of the syntax used by the main function: leave_COMPLETE:

function [sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_COMPLETE(y,id,firmid,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,filename)

%% DESCRIPTION
%This function computes leave out estimates in two-way fixed effects 
%model and conducts inference as described in KSS.

%The mandatory input is person-year dataset that has to be sorted
%by workers' identifiers (id) and year (xtset id year in Stata). The function
%automatically performs computation of the largest connected set and leave
%out connected set. 

%This function is specifically written not only to replicate the results of
%KSS but also to be applied to any other two-way fixed effects model
%(student-teachers, patient-doctors, etc).

%Check webpage for updates and other comments. 
%% DESCRIPTION OF THE INPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.

%id: worker indicators. Dimensions: N* x 1

%firmid: firm indicators. Dimensions: N* x 1

%leave_out_level: string variable that takes three values:

%'obs': perform leave-out by leaving a person-year observation out (default)

%'matches': perform leave-out by leaving an entire person-firm match out.

%'workers': perform leave-out by leaving an entire worker's history out.

%Please note that 'matches' and 'workers' options are still in beta mode
%and need further testing.
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                    %---NON-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%controls:
%Matrix of controls with dimensions: N* x K. This matrix of controls must
%be appropriately defined by the user ex-ante. For instance, if the user 
%wants to include time effects then the user should include in the matrix 
%'controls' the set of dummy variables associated to a particular year
%effect, making sure to avoid potential collinearity issues. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%resid_controls: Binary. 

%If 1 and input "controls" is not empty, then the code partials out 
%the effects of these controls before performing leave out estimates.
%In particular, in Step 1, after performing AKM on the largest connected
%set, we partial out the effect of controls using the corresponding AKM
%estimates.

%Default is 0. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
%andrews_estimates: Binary.

% If 1, the code reports homoskedastic corrections as described of 
% Andrews (2008). Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
%eigen_diagno: Binary. 

% If 1, the code outputs the lindeberg condition and 
% eigenvalue ratio of theorem 1. The code will also output the 
% weak-id confidence intervals using the AM method described in the paper 
% Default is 0. The eigenvalue ratio is calculated via simulations.  

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%subsample_llr_fit: Binary. 
%If 1, the code performs the non-parametric fit needed to compute standard 
%errors in the strong identification case using a subsample of the data. 
%See the end of Appendix B for further details.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%restrict_movers: Binary. 
%If 1, the code performs leave-out estimation just by focusing on sample
%of movers.
% Default is 0.
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%filename: string. 
%Used to name the saved outputs.
%Default is 'leave_one_out_estimates';
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

%% DESCRIPTION OF THE OUTPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%sigma2_psi: leave-out variance of firm effects. 
%sigma_psi_alpha: leave-out covariance of firm, person effects. 
%sigma2_alpha: leave-out variance of person effects. 

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
