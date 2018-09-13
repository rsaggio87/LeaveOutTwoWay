function [sigma2_psi,sigma_psi_alpha,sigma2_alpha,SE_sigma2_psi,SE_sigma_psi_alpha,SE_sigma2_alpha] = leave_out_COMPLETE(y,id,firmid,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename)
%% Author: Raffaele Saggio
%Email: rsaggio@princeton.edu

%% DESCRIPTION
%This function computes leave out estimates in two-way fixed effects 
%model and conducts inference as described in KSS.

%The mandatory input is person-year dataset that has to be sorted
%by workers' identifiers (id) and year (xtset id year in Stata). The function
%automatically performs computation of the largest connected set and leave
%out connected set. 

%This function can be applied to any two-way fixed effects model
%(student-teachers, patient-doctors, etc) and supports unbalanced data as
%well as T>2 panels.

%We use AKM jargon (workers, firms) when describing the code for
%simplicity.
%% Version:
% 1.0: Wrote documentation. 06.15.2018.

% 1.1:  Speed up computation of eigenvalues/vectors - by avoiding storage of
%       large matrices. 06.18.2018.

% 1.11: Eliminated the option Ndiagno. 06.19.2018.

% 1.12: Read Nargout and compute only parameters asked by user and added
%       option for whether user wants computation of standard error.
%       06.22.2018

% 1.2:  Added options to run fast Local Linear Regression, helpful in large
        %datasets. 01.07.2018 

% 1.25: Fix some bugs arising when setting leave_out_level to either
%       "matches" or "workers"  10.07.2018 

% 1.3:  Dropped stayers that have only one person year observations (for
%       which Pii=1 when estimating the model in levels). 25.07.2018

% 1.32: Improved readability of saved results at the end of the code. 
%       Changed locations of saved results 
%       Added export to .csv results 31.07.2018

% 1.5:  Significantly improved computation of Leave out matrices when there
%       are no controls in the model or the controls have been residualized
%       in a prior step. 
%          
%               We introduced the following changes:
%
%                1. Added CMG routine to speed computation of linear system
%                involving the Laplacian matrix as design matrix. Users need
%                to make sure that the CMG package is installed and placed
%                in the main directory as shown in the GitHub repository.
%
%                2. Read movers-stayers structure to fasten computation of (Bii,Pii).
%
%                In terms of speed, for the test dataset used in "example.m": 
%
%                1. With version 1.32 the code takes 260 seconds to compute (Bii,Pii).
%                2. With version 1.5 the code takes  23 seconds to compute (Bii,Pii).
% 
% 
% 1.51: Better management of large sparse matrices when invoking parfor to 
%       compute using the option parallel.pool.Constant.
%
% 1.52: Added more outputs to the function to simplify possible post-estimation 
%       commands.
%
% 1.55: Added more options to run non-parametric fit.
%
% 2.0:  Added option to approximate (Bii,Pii) using Random Projections methods
%       that build on the Johnson Lindestrauss Lemma - See Appendix B of KSS.
%
%       This especially helpful in massive datasets where exact computation 
%       of (Bii,Pii), even after the improvements introduced from version 1.5, 
%       is close to be prohibitive in terms of computation time.
%          
%       In particular, with this new release we added the following inputs:
%
%                1. type_of_algorithm: This takes two values: "exact" or "JLL".
%
%                    "exact": corresponds to exact computation of (Bii,Pii)
%                     as in version 1.5.          
%                    
%                    "JLL": applies random projection methods to
%                    approximate (Bii,Pii) as detailed in Appendix B of
%                    KSS. 
%                   
%                    Default is "exact".
%                
%                 2. epsilon: this governs the tradeoff b/w speed and unbiasdness 
%                    when estimating (Bii,Pii). Smaller values of epsilon implies 
%                    less bias but more computation time.
%
% 2.05: Made additional improvements for computation of (Bii,Pii) in "eff_res" by looking
%		at unique matches.
%
%       Note: The algorithm introduced with from release 2.0 only applies to
%       the case when there are no controls in the model or the controls 
%       have been residualized in a prior step so that the design matrix 
%       corresponds to a Laplacian matrix.     
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
%% 					DESCRIPTION OF THE INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.
%--
%id: worker indicators. Dimensions: N* x 1
%--
%firmid: firm indicators. Dimensions: N* x 1
%--
%leave_out_level: string variable that takes two values:

    %'obs': perform leave-out by leaving a person-year observation out (default)

    %'matches': perform leave-out by leaving an entire person-firm match out.
    
%Option 'matches' is currently in beta mode and needs further testing.

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
%resid_controls: 0 1 2. 

%If 0, the code includes the vector of controls specified by the user
%(provided that is not empty) in estimation of the two-way model.

%If 1 and input "controls" is not empty, then the code partials out 
%the effects of these controls before performing leave out estimates.
%In particular, in Step 1, after performing AKM on the largest connected
%set, we partial out the effect of controls using the corresponding AKM
%estimates.

%If 2 and input "controls" is not empty, the code acts as if the controls
%are not used in estimation of the model but saves the corresponding vector
%of controls in the leave out connected set.

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
% (assuming q=1) provided that do_SE=1.

% Default is 0.  

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%subsample_llr_fit: Discrete. 

%This governs how to compute the local linear regressions (LLR) fit needed
%to derive the standard errors as described in KSS when q=0. 
%See Remark 10 and appendix B of KSS.

%If 0, the code computes LLR (i.e. Lowess model) using the entire data 
%(i.e. both movers and stayers). 

%If 1, and there are no controls in the model, the code distinguishes
%between movers and stayers. For stayers, the assigned fitted value 
%is the average of \hat{\sigma}_i within cells defined in
%terms of unique values of T_i, where T_i is the number of person-year
%observations in which we observe a given worker. For movers, the code runs
%a standard LLR model. 

%If 1, and there are controls in the model, the code runs a stratified
%LLR separately for movers and stayers.

%If 2, and there are no controls in the model, the code creates "Kgrid"
%equally sized bins of Bii and Pii for movers only. 
%The code then fits a weighted LLR model across these "KGrid" x "KGrid" cells 
%weighting by cell size. For stayers, the assigned fitted value 
%is the average of \hat{\sigma}_i within cells defined in
%terms of unique values of T_i, where T_i is the number of person-year
%observations in which we observe a given worker. 

%If 2, and there are controls in the model, the code creates "Kgrid"
%equally sized bins of Bii and Pii for both movers and stayers. 
%We then fit a stratified binned LLR model separately 
%for movers and stayers across these "KGrid" x "KGrid" cells weighting 
%the estimates by cell size. 

%If 3, the code creates "Kgrid" equally sized bins of Bii and Pii
%The code then fits a weighted LLR model across these "KGrid" x "KGrid" 
%cells weighting by cell size. 

%If 4, the code fits a multivariate Nadaraya Watson style of non parametric
%estimator by averaging the value of \hat{\sigma}_i within "KGrid" x "KGrid" 
%cells of (Bii,Pii).

%See the function "llr_fit" for further details. Default value of
%"Kgrid=1000".

%Default value is 2. 

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%restrict_movers: Binary. 
%If 1, the code performs leave-out estimation just by focusing on sample
%of movers.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%do_SE: do_SE. 
%If 1, the code provides calculations of the standard errors assuming
%strong identification (q=1). 
% Default is 1.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%type_of_algorithm: 
%This takes two values: "exact" or "JLL".                    

%   "exact": performs exact computation of (Bii,Pii). 

%   "JLL": perform random projection methods to approximate (Bii,Pii) as 
%   detailed in Appendix B of KSS. 

%In large datasets (say where #of Workers + #of Firms > 50K) we suggest
%setting type_of_algorithm="JLL". As described in Table B1 in large
%datasets JLL gives back idenentical estimates compared to exact methods
%while taking 20-30 times less computation time. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%epsilon: This should be b/w 0 and 1.

%epsilon governs the tradeoff b/w speed and accuracy when estimating 
%(Bii,Pii) using type_of_algorithm=JLL.

%Smaller values of epsilon implies less accuracy but higher computation time.

%Default is 0.1. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%filename: string. 
%Where saved results should be stored and named. Use name like 
%"leave_out_results" and not "leave_out_results.csv"

%Default is 'leave_one_out_estimates';
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

% DESCRIPTION OF THE OUTPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%sigma2_psi: leave-out variance of firm effects.
%SE_sigma2_psi: Standard Error (SE) of leave-out variance of firm effects
%(q=0, high rank case)

%sigma_psi_alpha: leave-out covariance of firm, person effects.
%SE_sigma_psi_alpha: SE of leave-out covariance of firm, person effects.

%sigma2_alpha: leave-out variance of person effects. 
%SE_sigma2_alpha: SE of leave-out variance of person effects.

%y: outcome variable in the leave out connected set.

%id: re-normalized worker identifiers in the leave out connected set.

%firmid: re-normalized firm identifiers in the leave out connected set.

%controls: vector of extra controls in the leave out connected set.

%Pii: (diagonal) statistical leverage associated with two way model.

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

% SAVED FILES

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%The code saves the following files:

%One .mat file after step 3  with all matrices in memory.
%One .mat file after completing the code with all matrices in memory.
%One .mat file containing the relevant leave-out leverages.
%A   .csv file with location and name specified by the user.

%The .csv saves the following variables belonging 
%to the leave out connected set.

%y: outcome variable.
%firmid: firm identifier. (normalized)
%id: worker identifier. (normalized)
%firmid_old: firm identifier. (as in original input data)
%firmid: firm identifier. (as in original input data)
%controls: vector of controls
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%% READ
no_controls=0;

if nargin < 4
    error('More arguments needed');
end

if nargin == 4 
    no_controls=1;
    controls=ones(size(y,1),1);
    resid_controls=0;
    andrews_estimates=0;
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.1;   
    filename='leave_one_out_estimates';
end

if nargin == 5   
    no_controls=1;
    resid_controls=1; 
    controls=ones(size(y,1),1);
    andrews_estimates=0;
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;   
    filename='leave_one_out_estimates';
end

if nargin == 6
    andrews_estimates=0;    
    eigen_diagno=0;
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 7    
    eigen_diagno=0;    
    subsample_llr_fit=0;
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 8      
    subsample_llr_fit=0;    
    restrict_movers=0;
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 9     
    restrict_movers=0; 
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 10  
    do_SE=1;
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 11 
    type_of_algorithm='exact';
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 12 
    epsilon=0.005;
    filename='leave_one_out_estimates';
end

if nargin == 13 
    filename='leave_one_out_estimates';
end


if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
    resid_controls=0;
end

if size(controls,2)>0 && resid_controls==2
    no_controls=1;
    resid_controls=0;
end


if resid_controls==1 && no_controls== 1 
    error('cannot residualize if there are no controls specified')    
end

if resid_controls==0 && no_controls== 0 && strcmp(type_of_algorithm,'JLL')
    error('cannot run JLL algorithm on non-Laplacian design matrix. Please set resid_controls=1')    
end


%Read number of outputs
if  nargout==1
    n_of_parameters=1;
end

if nargout==2
    n_of_parameters=2;
end

if nargout>=3
    n_of_parameters=3;
end


%Read possible inconsistencies
if ~strcmp(leave_out_level,'obs') && restrict_movers == 0
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
disp('WARNING!!!!!!!!!!!!!!!')
disp('Variance of Person Effects not identified for stayers when leaving more than a single observation out.')
disp('The code will omit all output concerning the variance of person effects.')
disp('Consider setting "restrict_movers=1" if user is interested in the variance of person effects.')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp(s)
n_of_parameters=2;
end

%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
n_of_parameters
leave_out_level
no_controls
resid_controls
andrews_estimates
eigen_diagno
subsample_llr_fit
restrict_movers
do_SE
type_of_algorithm
epsilon
filename
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)





%% STEP 1: PRELIMINARIES
%As first step in our analysis, we run estimation of a standard AKM model
%on the original input data. 

%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Find the connected set. Then define dimensions and relevant matrices.
[y,id,firmid,id_old,firmid_old,controls] = connected_set(y,id,firmid,lagfirmid,controls);


%Define
NT=size(y,1);
J=max(firmid);
N=max(id);
D=sparse(1:NT,id',1);
F=sparse(1:NT,firmid',1);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
if no_controls==1
X=[D,F*S];
end
if no_controls==0
X=[D,F*S,controls];
end

%Run AKM
disp('Running AKM...')
tic
xx=X'*X;
xy=X'*y;
tic
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
toc
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
toc
ahat=b(1:N);
ghat=b(N+1:N+J-1);
pe=D*ahat;
fe=F*S*ghat;
xb=X*b;
r=y-xb;
dof=NT-size(X,2)-1;
TSS=sum((y-mean(y)).^2);
R2=1-sum(r.^2)/TSS;
adjR2=1-sum(r.^2)/TSS*(NT-1)/dof;


%Some auxiliaries for summary statistics.
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker
stayer=(firmid==lagfirmid);
stayer(gcs==1)=1;
stayer=accumarray(id,stayer);
T=accumarray(id,1);
stayer=T==stayer;
movers=stayer~=1;
movers=movers(id);
id_movers=id(movers);
[~,~,n]=unique(id_movers);
Nmovers=max(n);

%Run AKM
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['STEP 1: AKM Estimates on Largest Connected Set '];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(max(Nmovers))];
disp(s);
s=['# of Firms: ' num2str(size(F,2))];
disp(s);
s=['# of Person Year Observations: ' num2str(sum(diag((F'*F))))];
disp(s);
s=['-*-*-*-*-*-*AKM RESULTS-*-*-*-*-*-*'];
disp(s)
COV=cov(fe,pe);
s=['Variance of Firm Effects: ' num2str(COV(1,1))];
disp(s);
s=['Covariance of Firm and Person Effects: ' num2str(COV(1,2))];
disp(s);
s=['Variance of Person Effects: ' num2str((COV(2,2)))];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str((corr(fe,pe)))];
disp(s);
s=['R2: ' num2str(R2)];
disp(s);
s=['Adj.R2: ' num2str(adjR2)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%Do some cleaning of matrices in memory
clear xx xy L xb pe fe ahat ghat F D S Lchol


%%Deresidualized outcome variable (allows the leave one out methodology to run faster)
if resid_controls==1
y=y-X(:,N+J:end)*b(N+J:end);
no_controls=1;
end
clear X b r

%Call CMG routine if there are no controls from the model
if no_controls == 1
warning('off', 'all') 
cd codes;    
path(path,'CMG'); %this contains the main LeaveOut Routines.
MakeCMG;
cd ..
end

%% STEP 2: LEAVE ONE OUT CONNECTED SET
%Here we compute the leave out connected set as defined in Appendix B. 
%The input data is represented by the largest connected set. After applying
%the function 'pruning_unbal_v3', the output data will be a connected set
%such that the associated bipartite graph between workers and firms remains
%connected after removing any particular worker from the graph.


%Focus estimation on movers only?
if restrict_movers==1
sel=movers;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%reset ids
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n;
end

%Leave One Out Largest Connected Set
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[y,firmid,id,id_old,firmid_old,controls] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%%%Drop stayers with a single person year observation
T=accumarray(id,1);
T=T(id);
sel=T>1;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);

%Resetting ids one last time.
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n; 

%Important Auxiliaries
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker
stayer=(firmid==lagfirmid);
stayer(gcs==1)=1;
stayer=accumarray(id,stayer);
T=accumarray(id,1);
stayer=T==stayer;
movers=stayer~=1;
movers=movers(id);
T=T(id);
id_movers=id(movers);
[~,~,n]=unique(id_movers);
Nmovers=max(n);
NT=size(y,1);
D=sparse(1:NT,id',1);
N=size(D,2);
F=sparse(1:NT,firmid',1);
J=size(F,2);
if no_controls==0
K=size(controls,2);
end
if no_controls==1
K=0;
end
clear id_mover
%Summarize
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(max(firmid))];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)

%% STEP 3: COMPUTATION OF (Pii,Bii)
%The code proceeds as follows: 

%STEP 3A:
%We start by reading at what level the user wants to perform leave out 
%(person-year, match, person). Based on this we create a list of 
%observations that has to be jointly removed, see the function 
%'index_constr'.

%STEP 3B:
%We then build the design matrix as a grounded Laplacian matrix.

%STEP 3C:
%Next we compute estimation of the leave out matrices for each variance
%decomposition parameter, using the function 'eff_res'. 


%% STEP 3A: Indexes needed to construct Bii, Pii
[~,~,match_id]=unique([id firmid],'rows'); 
if strcmp(leave_out_level,'obs')
index=(1:NT)';    
clustering_var=index;
end
if strcmp(leave_out_level,'matches')
clustering_var=match_id;
end
elist=index_constr(clustering_var,id,match_id);

%% STEP 3B: Build Design
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];

%Design matrix
if no_controls==0
X=[D,F*S,controls];
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
end
if no_controls==1
X=[D,-F];
xx=X'*X;
disp('Building preconditioner for Laplacian Matrix...')
Lchol = cmg_sdd(xx); %preconditioner for Laplacian matrices.
end

%% STEP 3C: NOW COMPUTE THE LEAVE OUT MATRICES
tic
disp('Calculating Leave Out Matrices...')

if n_of_parameters==1
    [Lambda_P, Lambda_B_fe]=eff_res(X,xx,Lchol,N,J,K,elist,leave_out_level,movers,T,type_of_algorithm,id,firmid,epsilon);
end

if n_of_parameters==2
    [Lambda_P, Lambda_B_fe,Lambda_B_cov]=eff_res(X,xx,Lchol,N,J,K,elist,leave_out_level,movers,T,type_of_algorithm,id,firmid,epsilon);
end

if n_of_parameters==3
    [Lambda_P, Lambda_B_fe,Lambda_B_cov,Lambda_B_pe]=eff_res(X,xx,Lchol,N,J,K,elist,leave_out_level,movers,T,type_of_algorithm,id,firmid,epsilon);
end

clear elist index clustering_var

disp('Time to compute Leave one out matrices')
toc
disp('Saving (Bii,Pii)')
s=[filename '_after_step3'];
save(s)

%Reshape the X back to grounded Laplacian.
if K == 0
    X=[D,F*S];
    xx=X'*X;
    Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
end

%% STEP 4: DIAGNOSTICS
%This part is only computed provided that the option 'eigen_diagno' is
%turned on. In this part, we calculate the squared eigenvalue ratio and
%Lindeberg conditions that constitute the key conditions to verify the
%validity of Theorem 1. Lindeberg condition assumes q=1.

if eigen_diagno==1 
EIG_NORM=zeros(3,n_of_parameters);
max_x1bar_sq=zeros(n_of_parameters,1);
lambda_1=zeros(n_of_parameters,1);
SUM_EIG=zeros(n_of_parameters,1);
x1bar_all=zeros(NT,n_of_parameters);
tic    
%Begin by calcuting the sum of the squared of the eigenvalues for the corresponding Atilde
disp('Calculating Sum of Squared Eigenvalues via Hutchinson...')

if n_of_parameters==1
    [trace_fe]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==2
    [trace_fe, trace_cov]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==3
    [trace_fe,trace_cov,trace_pe]=trace_Atilde_sqr(X,F*S,D,xx,Lchol);
end

if n_of_parameters==1
    SUM_EIG(1)=trace_fe;
end

if n_of_parameters==2
    SUM_EIG(1)=trace_fe;
    SUM_EIG(2)=trace_pe;
end

if n_of_parameters==3
    SUM_EIG(1)=trace_fe;
    SUM_EIG(2)=trace_pe;
    SUM_EIG(3)=trace_cov;
end
disp('done') 

%Issue: we use a trick in order to compute the eigenvalues/vectors
%associated with the matrix Atilde without having to store the matrix
%Atilde in memory

for pp=1:n_of_parameters
  
    if pp == 1
    type_quadratic_form='fe';
    entry=['Calculating Diagnostic for Variance of Firm Effects...'];
    end
    
    if pp == 2
    type_quadratic_form='cov';
    entry=['Calculating Diagnostic for Covariance of Firm, Person Effects...'];
    end
    
    if pp == 3
    type_quadratic_form='pe';
    entry=['Calculating Diagnostic for Variance of Person Effects...'];
    end
    
    %Calculate
    disp(entry)
    [Q,lambda_eig] = eigAux(type_quadratic_form,xx,Lchol,F*S,D,K);
    lambda_eig=diag(lambda_eig);
    lambda_1(pp)=lambda_eig(1); 
    [EIG_NORM(:,pp),x1bar_all(:,pp)] = eig_x1bar(X,Q,lambda_eig,SUM_EIG(pp)); %Note: What we label as x1bar in the code corresponds to w_i1 in KSS.
    max_x1bar_sq(pp)=max(x1bar_all(:,pp).^2);
    
    
end
disp('checking, must report 1')
sum(x1bar_all.^2,1)  
end

%% STEP 5: ESTIMATION
%We now conduct estimation on leave-out connected set. 

%These are the steps:

%Step 5.1: Here we perform standard AKM (i.e. plug-in) estimates of the
%variance decomposition parameters.


%Step 5.2: Here we perform leave-out estimates and inference of the
%variance decomposition parameters.

%Step 5.3: Here we perform homoskedasticity corrected estimates of the
%variance decomposition parameters.



%% STEP 5.1: AKM (PLUG-IN)
xy=X'*y;
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
ahat=b(1:N);
ghat=b(N+1:N+J-1);
pe=D*ahat;
fe=F*S*ghat;
xb=X*b;
eta=y-xb;
dof=NT-size(X,2)-1;
TSS=sum((y-mean(y)).^2);
R2=1-sum(eta.^2)/TSS;
adjR2=1-sum(eta.^2)/TSS*(NT-1)/dof;
COV=cov(fe,pe);
sigma_2_psi_AKM=COV(1,1);
sigma_2_alpha_AKM=COV(2,2);
sigma_alpha_psi_AKM=COV(1,2);

%% STEP 5.2: Leave One Out Estimation.
I_Lambda_P=(speye(NT,NT)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
msg = lastwarn ; 
if ~isempty(strfind(msg, 'singular')) 
            	s=['******************************************'];
                disp(s);
                disp(s); 
				disp('Warning: OLS coefficient not always identified when leaving a particular set of observation out as specified by "leave_out_level"')
				disp('One example where this occurs is when the user asks to run leave out on matches without restricting the analysis to movers only.')
				s=['******************************************'];
                disp(s);
                disp(s); 
				
end
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));


%Preallocate vectors with results
theta=zeros(n_of_parameters,1);
V_theta=zeros(n_of_parameters,1);
COV_R1=zeros(2,2,n_of_parameters);
gamma_sq=zeros(n_of_parameters,1);
F_stat=zeros(n_of_parameters,1);
b_1=zeros(n_of_parameters,1);
theta_1=zeros(n_of_parameters,1);

%Loop over parameters to be estimated
for pp=1:n_of_parameters
    
    if eigen_diagno == 1
        x1_bar=x1bar_all(:,pp);
    end
    
    if pp == 1 %Variance of Firm Effects
        type_quadratic_form='fe';
        Lambda_B=Lambda_B_fe;
        bias_part=sigma_2_psi_AKM; 
        if no_controls == 0
            A_b=[zeros(N,1); S'*F'*(fe-mean(fe)); zeros(K,1)];    
        end 
        if no_controls == 1
            A_b=[zeros(N,1); S'*F'*(fe-mean(fe))];
        end
        entry=['Calculating Leave out, Variance of Firm Effects...'];
    end


    if pp == 2 %CoVariance of Person,Firm Effects
        type_quadratic_form='cov';
        Lambda_B=Lambda_B_cov;
        bias_part=sigma_alpha_psi_AKM; 
        if no_controls==0
        A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe)); zeros(K,1)];
        end
        if no_controls==1
        A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe))];
        end
        entry=['Calculating Leave out, Covariance of Person, Firm Effects...'];
    end

    if pp == 3 %Variance of Person Effects
        type_quadratic_form='pe';
        Lambda_B=Lambda_B_pe;
        bias_part=sigma_2_alpha_AKM; 
        if no_controls==0
            A_b=[D'*(pe-mean(pe)); zeros(J-1,1); zeros(K,1)];
        end
        if no_controls==1
            A_b=[D'*(pe-mean(pe)); zeros(J-1,1)];
        end
        entry=['Calculating Leave out, Variance of Person Effects...'];
    end

%Tell me what are we doing
     disp(entry)
   
%Auxiliary
    if do_SE  == 1 
        [W_to_use, my_first_part] = construc_W(y,X,xx,Lchol,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);
    end
%Compute the non-parametric fit needed for standard error estimation in 
%high rank case.
    if do_SE == 1  
        sigma_predict= llr_fit(Lambda_P,Lambda_B,y,eta_h,subsample_llr_fit,K,movers,T);
    end

%Run
    if eigen_diagno == 1  && do_SE == 0
        [theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B);
    end

    if eigen_diagno == 0  && do_SE == 0
        [theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B);
    end
    
    if eigen_diagno == 0  && do_SE == 1
        [theta(pp), V_theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict);
    end
    if eigen_diagno == 1  && do_SE == 1
        [theta(pp), V_theta(pp), COV_R1(:,:,pp), gamma_sq(pp), F_stat(pp), b_1(pp), theta_1(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict,x1_bar,lambda_1(pp));
    end
    
    

end

clear A_b Lambda_B W_to_use I_Lambda_P L_P

if n_of_parameters==1 && do_SE==0
    sigma2_psi=theta(1);
end

if n_of_parameters==2 && do_SE==0
    sigma2_psi=theta(1);
    sigma_psi_alpha=theta(2);
end   

if n_of_parameters==3 && do_SE==0
    sigma2_psi=theta(1);
    sigma_psi_alpha=theta(2);
    sigma2_alpha=theta(3);
end


if n_of_parameters==1 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
end

if n_of_parameters==2 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
    sigma_psi_alpha=theta(2);
    SE_sigma_psi_alpha=sqrt(V_theta(2));
end   

if n_of_parameters==3 && do_SE==1
    sigma2_psi=theta(1);
    SE_sigma2_psi=sqrt(V_theta(1));
    sigma_psi_alpha=theta(2);
    SE_sigma_psi_alpha=sqrt(V_theta(2));
    sigma2_alpha=theta(3);
    SE_sigma2_alpha=sqrt(V_theta(3));
end  

%R2
if n_of_parameters==3
    explained_var_leave_out=(sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var(y);
end




%% %% STEP 5.3: Andrews
if andrews_estimates==1 && no_controls==0    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D,controls);
end
if andrews_estimates==1 && no_controls==1    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D);
end


%% STEP 6: REPORTING
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
s=['Results: Leave One Out Largest Connected Set'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(size(F,2))];
disp(s);
s=['# of Person Year Observations: ' num2str(NT)];
disp(s);
s=['Maximum Leverage: ' num2str(max(diag(Lambda_P)))];
disp(s);
s=['-*-*-*-*-*-*AKM*-*-*-*'];
disp(s)
s=['Variance of Firm Effects: ' num2str(sigma_2_psi_AKM)];
disp(s)
s=['Covariance of Firm and Person Effects: ' num2str(sigma_alpha_psi_AKM)];
disp(s);
s=['Variance of Person Effects: ' num2str(sigma_2_alpha_AKM)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str((corr(fe,pe)))];
disp(s);
s=['R2: ' num2str(R2)];
disp(s);
if andrews_estimates == 1
s=['-*-*-*-*-*-*ANDREWS*-*-*-*'];
disp(s);
s=['Variance of Firm Effects: ' num2str(var_corrected_fe)];
disp(s);
s=['Covariance of Firm and Person Effects: ' num2str(var_corrected_cov)];
disp(s);
s=['Variance of Person Effects: ' num2str(var_corrected_pe)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str(var_corrected_cov/(sqrt(var_corrected_fe)*sqrt(var_corrected_pe)))];
disp(s);
s=['Total Explained Variation: ' num2str(adjR2)];
disp(s);
end
s=['-*-*-*-*-*-*LEAVE ONE OUT*-*-*-*'];
disp(s)
s=['Variance of Firm Effects: ' num2str(sigma2_psi)];
disp(s);
if n_of_parameters>=2
s=['Covariance of Firm and Person Effects: ' num2str(sigma_psi_alpha)];
disp(s);
end
if n_of_parameters==3
s=['Variance of Person Effects: ' num2str(sigma2_alpha)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
disp(s);
s=['Total Explained Variation: ' num2str(explained_var_leave_out)];
disp(s);
end
%% Alternative Way to Show Output
s=['-*-*-*-*-*-*Variance of Firm Effects-*-*-*-*-*-*'];
disp(s)
s=['AKM: ' num2str(sigma_2_psi_AKM)];
disp(s)
if andrews_estimates == 1
s=['Andrews: ' num2str(var_corrected_fe)];
disp(s)
end
s=['Leave One Out: ' num2str(sigma2_psi)];
disp(s)
if do_SE == 1
s=['Leave One Out SE: ' num2str(SE_sigma2_psi)];
disp(s)
end
if n_of_parameters>=2
    s=['-*-*-*-*-*-*Covariance of Firm,Person Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str(sigma_alpha_psi_AKM)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_cov)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma_psi_alpha)];
    disp(s);
    if do_SE==1
    s=['Leave One Out SE: ' num2str(SE_sigma_psi_alpha)];
    disp(s)
    end
end
if n_of_parameters>=3
    s=['-*-*-*-*-*-*Variance of Person Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str(sigma_2_alpha_AKM)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_pe)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma2_alpha)];
    disp(s);
    if do_SE==1
    s=['Leave One Out SE: ' num2str(SE_sigma2_alpha)];
    disp(s)
    end
    s=['-*-*-*-*-*-*Correlation of Person,Firm Effects-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: ' num2str((corr(fe,pe)))];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(var_corrected_cov/(sqrt(var_corrected_fe)*sqrt(var_corrected_pe)))];
    disp(s);
    end
    s=['Leave One Out: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
    disp(s);
    s=['-*-*-*-*-*-*R2-*-*-*-*-*-*'];
    disp(s)
    s=['AKM: '  num2str(R2)];
    disp(s);
    if andrews_estimates == 1
    s=['Andrews: ' num2str(adjR2)];
    disp(s);
    end
    s=['Leave One Out: ' num2str(explained_var_leave_out)];
    disp(s);
end
%% Focus on Diagnostics
%Note: The eigenvalue ratios (and sum of squared eigenvalues) 
%      reported in Table 3 of KSS have been calculated exactly whereas here
%      we report the results obtained via simulations. 
%      The differences are neglible as one can see from below.
if eigen_diagno==1
    for pp=1:n_of_parameters
        
        if pp == 1
            title='Diagnostics on Variance of Firm Effects';
        end  
        if pp == 2
            title='Diagnostics on CoVariance of Person, Firm Effects';
        end
        if pp == 3
            title='Diagnostics on Variance of Person Effects';
        end

        s=['*********************' title '*********************'];
        disp(s);
        s=['ratio of eigenvalues: '];
        disp(s)
        EIG_NORM(1:3,pp)
        s=['Weak Lindeberg Condition: ' num2str(max_x1bar_sq(pp))];
        disp(s)
        s=['Sum of squared eigenvalues: ' num2str(SUM_EIG(pp)/NT^2)];
        disp(s)
        s=['Variance of b_1:  ' num2str(COV_R1(1,1,pp))];
        disp(s);
        s=['Variance of \hat{\theta}_1: ' num2str(COV_R1(2,2,pp))];
        disp(s);
        s=['Covariance of (b_1,\theta_1): ' num2str(COV_R1(1,2,pp))];
        disp(s);
        s=['Correlation of (b_1,\theta_1): ' num2str((COV_R1(1,2,pp))/(sqrt(COV_R1(2,2,pp))*sqrt(COV_R1(1,1,pp))))];
        disp(s);
        s=['gamma squared: ' num2str(gamma_sq(pp))];
        disp(s);
        s=['Fstatistic: ' num2str(F_stat(pp))];
        disp(s);
        s=['******************************************'];
        disp(s);
    end
end
%% Focus on Inference
if do_SE==1
    for pp=1:n_of_parameters
        
        if pp == 1
            title='Inference on Variance of Firm Effects';
        end
        
        if pp == 2
            title='Inference on CoVariance of Person, Firm Effects';
        end
        
        if pp == 3
            title='Inference on Variance of Person Effects';
        end
        
        s=['*********************' title '*********************'];
        disp(s);
        s=['SE under q=0: ' num2str(sqrt(V_theta(pp)))];
        disp(s);
        s=['CI under q=0: ' num2str(theta(pp)-1.96*sqrt(V_theta(pp))) '  '  '  ' num2str(theta(pp)+1.96*sqrt(V_theta(pp)))];
        disp(s);   
        
        if eigen_diagno==1
                [UB,LB,C]=AM_CI(NT,lambda_1(pp),gamma_sq(pp),COV_R1(:,:,pp),b_1(pp),theta_1(pp));
                s=['CI under q=1: ' num2str(LB) '  '  '  ' num2str(UB)];
                disp(s);  
                s=['Curvature: ' num2str(C)];
                disp(s); 
                s=['******************************************'];
                disp(s); 
        end
    end
end
%Save Outputs.
s=[filename '_completed'];
save(s)
s=[filename '_Lambda_P'];
save(s,'Lambda_P')


%Export csv with data from leave out connected set
Pii=full(diag(Lambda_P));
out=[y,firmid,id,firmid_old,id_old,Pii,controls];
s=[filename '.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
end

