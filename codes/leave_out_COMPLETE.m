function [sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_COMPLETE(y,id,firmid,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,Ndiagno,filename)
%% Author: Raffaele Saggio
%Email: raffaele.saggio@berkeley.edu

%% Version:
% 1.0: Wrote documentation. 06.15.2018.


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
%Ndiagno: Takes three values. 
%1: Code reports conditions of Theorem 1 for Variance of Firm Effects only. 
%2: Code reports conditions of Theorem 1 for Variance of Firm Effects and 
    %Covariance of Person and Firm effects. 
%3: Code reports conditions of Theorem 1 for Variance of Firm Effects, 
%   %Covariance of Person and Firm effects and Variance of person effects.

%Default is 1.
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

%% SAVED FILES
%The code saves the following files:

%One mat file after step 3.
%One mat file after finishing computation with all the relevant matrices.
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
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 5
resid_controls=1;    
andrews_estimates=0;
eigen_diagno=0;
subsample_llr_fit=0;
restrict_movers=0;
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 6
andrews_estimates=0;    
eigen_diagno=0;
subsample_llr_fit=0;
restrict_movers=0;
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 7    
eigen_diagno=0;    
subsample_llr_fit=0;
restrict_movers=0;
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 8      
subsample_llr_fit=0;    
restrict_movers=0;
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 9     
restrict_movers=0; 
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 10  
Ndiagno=1;
filename='leave_one_out_estimates';
end

if nargin == 11  
Ndiagno=1;
filename='leave_one_out_estimates';
end


if size(controls,2)==0
no_controls=1;
controls=ones(size(y,1),1);
resid_controls=0;
end


if resid_controls==1 && no_controls== 1 
error('cannot residualize if there are no controls specified')    
end
%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
leave_out_level
no_controls
resid_controls
andrews_estimates
eigen_diagno
subsample_llr_fit
restrict_movers
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
clear X
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
%Do some cleaning of matrices in memory
clear xx xy b L xb pe fe ahat ghat F D S Lchol


%%Deresidualized outcome variable (allows the leave one out methodology to run faster)
if resid_controls==1
y=y-X(:,N+J:end)*b(N+J:end);
no_controls=1;
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

%%%Resort data if running leave one out on matches
if strcmp(leave_out_level,'matches')
[~,~,match_id]=unique([id firmid],'rows');  
[~,IX]=sort(match_id);
y=y(IX);
firmid=firmid(IX);
id=id(IX);
id_old=id_old(IX);
firmid_old=firmid_old(IX);
controls=controls(IX,:);
match_id=match_id(IX);
end


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
if strcmp(leave_out_level,'workers') 
clustering_var=id;
end
if strcmp(leave_out_level,'obs')
index=(1:NT)';    
clustering_var=index;
end
if strcmp(leave_out_level,'matches')
clustering_var=match_id;
end
elist=index_constr(clustering_var,id);

disp('checking, must report 1')
mean(elist(:,1)<=elist(:,2))

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
X=[D,F*S];
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
end

%% STEP 3C: NOW COMPUTE THE LEAVE OUT MATRICES
tic
[Lambda_P, Lambda_B_fe, Lambda_B_cov, Lambda_B_pe]=eff_res(X,xx,Lchol,N,J,elist,leave_out_level);
disp('Time to compute Leave one out matrices')
toc

s=['mat/after_step_3_complete' filename];
save(s)


disp('checking, must report 1')
issymmetric(Lambda_P)
issymmetric(Lambda_B_fe)
issymmetric(Lambda_B_cov)
issymmetric(Lambda_B_pe)

%% STEP 4: DIAGNOSTICS
%This part is only computed provided that the option 'eigen_diagno' is
%turned on. In this part, we calculate the squared eigenvalue ratio and
%Lindeberg conditions that constitute the key conditions to verify the
%validity of Theorem 1. Lindeberg condition assumes q=1.

%Future releases of this function will allow for faster computation of the 
%eigenvalue conditions. In particular, the step where we build and store 
%the A matrix for calculation of the eigenvectors can be sidestepped.

if eigen_diagno==1 
x1bar_pe=zeros(NT,1);
x1bar_cov=zeros(NT,1);
EIG_NORM=zeros(3,3);
max_x1bar_sq=zeros(3,1);
lambda_1=zeros(3,1);
SUM_EIG=zeros(3,3);
tic    
%Keep track of #of py observations associated with each worker-firm
degree_workers=accumarray(id,1);
degree_firms=diag((F*S)'*(F*S));

%Begin by calcuting the sum of the squared of the eigenvalues for the corresponding Atilde
disp('Calculating Sum of Squared Eigenvalues via Hutchinson...')
[trace_fe,trace_pe,trace_cov] = trace_Atilde_sqr(X,F*S,D,xx,Lchol);
SUM_EIG(1)=trace_fe;
SUM_EIG(2)=trace_pe;
SUM_EIG(3)=trace_cov;
disp('done') 

%Variance of firm effects
disp('Calculating Diagnostic for Variance of Firm Effects...')
degree=(F*S)'*F*S-(1/NT)*(degree_firms)*degree_firms'; 
if K>0
degree=[sparse(N,N+J+K-1); sparse(J-1,N) degree sparse(J-1,K); sparse(K,K+N+J-1)];
end
if K==0
degree=[sparse(N,N+J-1); sparse(J-1,N) degree];
end 
[Q, lambda_eig] = eigs(degree,xx,3);
clear degree
lambda_eig=diag(lambda_eig);
lambda_1(1)=lambda_eig(1);
%Calculate normalized eigenvalues and corresponding x1bar. 
[lambda_fe,x1bar_fe] = eig_x1bar(X,Q,lambda_eig,trace_fe);
EIG_NORM(:,1)=lambda_fe;
max_x1bar_sq(1)=max(x1bar_fe.^2);
disp('done')
%Note: What we label as x1bar in the code corresponds to w_i1 in KSS.

if Ndiagno == 3
%Variance of person effects.
disp('Calculating Diagnostic for Variance of Person Effects...')
degree=D'*D-(1/NT)*(degree_workers)*degree_workers'; 
if K>0
degree=[degree sparse(N,K+J-1); sparse(J+K-1,N+J+K-1)];
end
if K==0
degree=[degree sparse(N,J-1); sparse(J-1,N+J-1)];
end
[Q, lambda_eig] = eigs(degree,xx,3);
clear degree
lambda_eig=diag(lambda_eig);
lambda_1(2)=lambda_eig(1);
%Calculate normalized eigenvalues and corresponding x1bar.
[lambda_pe,x1bar_pe] = eig_x1bar(X,Q,lambda_eig,trace_pe);
EIG_NORM(:,2)=lambda_pe;
max_x1bar_sq(2)=max(x1bar_pe.^2);
disp('done.')
end 

if Ndiagno == 2
%Covariance of person, firm effects.
disp('Calculating Diagnostic for CoVariance of Person,Firm Effects...')
degree=D'*F*S-(1/NT)*(degree_workers)*degree_firms'; 
if K>0
degree=[sparse(N,N) 0.5*degree sparse(N,K); 0.5*degree' sparse(J-1,J-1)  sparse(J,K); sparse(K,K+N+J-1)];
end
if K==0
degree=[sparse(N,N) 0.5*degree; 0.5*degree' sparse(J-1,J-1)];
end
[Q, lambda_eig] = eigs(degree,xx,3);
clear degree
lambda_eig=diag(lambda_eig);
lambda_1(3)=lambda_eig(1);
%Calculate normalized eigenvalues and corresponding x1bar.
[lambda_cov,x1bar_cov] = eig_x1bar(X,Q,lambda_eig,trace_cov);
EIG_NORM(:,3)=lambda_cov;
max_x1bar_sq(3)=max(x1bar_cov.^2);
disp('done.')
end
%Finish up
disp('Time to Calculate Eigenvalues Conditions')
toc

disp('checking, must report 1')
sum(x1bar_fe.^2)
if Ndiagno==3
sum(x1bar_pe.^2)
end
if Ndiagno==2
sum(x1bar_cov.^2)  
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
disp('Running AKM...')
tic
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
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[eta_h, flag]=pcg(I_Lambda_P,eta,1e-5,1000,L_P,L_P'); %leave-out residuals

%Preallocate vectors with results
theta=zeros(3,1);
V_theta=zeros(3,1);
COV_R1=zeros(2,2,3);
gamma_sq=zeros(3,1);
F_stat=zeros(3,1);
b_1=zeros(3,1);
theta_1=zeros(3,1);

%Loop over parameters to be estimated
for pp=1:3
    
if pp == 1 %Variance of Firm Effects
type_quadratic_form='fe';
Lambda_B=Lambda_B_fe;
bias_part=sigma_2_psi_AKM; 
if eigen_diagno == 1
x1_bar=x1bar_fe;
end
if no_controls == 0
A_b=[zeros(N,1); S'*F'*(fe-mean(fe)); zeros(K,1)];    
end 
if no_controls == 1
A_b=[zeros(N,1); S'*F'*(fe-mean(fe))];
end
end

if pp == 2 %Variance of Person Effects
type_quadratic_form='pe';
Lambda_B=Lambda_B_pe;
bias_part=sigma_2_alpha_AKM; 
if eigen_diagno == 1
x1_bar=x1bar_pe;
end
if no_controls==0
A_b=[D'*(pe-mean(pe)); zeros(J-1,1); zeros(K,1)];
end
if no_controls==1
A_b=[D'*(pe-mean(pe)); zeros(J-1,1)];
end
end

if pp == 3 %CoVariance of Person,Firm Effects
type_quadratic_form='cov';
Lambda_B=Lambda_B_cov;
bias_part=sigma_alpha_psi_AKM; 
if eigen_diagno == 1
x1_bar=x1bar_cov;
end
if no_controls==0
A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe)); zeros(K,1)];
end
if no_controls==1
A_b=[0.5*D'*(fe-mean(fe)); 0.5*S'*F'*(pe-mean(pe))];
end
end

%Auxiliary
[W_to_use, my_first_part] = construc_W(y,X,xx,Lchol,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);

%Compute the non-parametric fit needed for standard error estimation in 
%high rank case.
sigma_predict= llr_fit(Lambda_P,Lambda_B,y,eta_h,subsample_llr_fit);

%Run
if eigen_diagno == 0 
[theta(pp), V_theta(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict);
end
if eigen_diagno == 1 
[theta(pp), V_theta(pp), COV_R1(:,:,pp), gamma_sq(pp), F_stat(pp), b_1(pp), theta_1(pp)]= leave_out_estimation_two_way(type_quadratic_form,y,J,N,X,xx,Lchol,bias_part,eta_h,I_Lambda_P,L_P,Lambda_B,W_to_use,my_first_part,sigma_predict,x1_bar,lambda_1(pp));
end
end
clear A_b Lambda_B W_to_use I_Lambda_P L_P
%Rename some things
sigma2_psi=theta(1);
sigma2_alpha=theta(2);
sigma_psi_alpha=theta(3);
V_sigma2_psi=V_theta(1);
V_sigma2_alpha=V_theta(2);
V_sigma_alpha_psi=V_theta(3);
%R2
explained_var_leave_out=(sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var(y);





%% %% STEP 5.3: Andrews
if andrews_estimates==1 && no_controls==0    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D,controls);
end
if andrews_estimates==1 && no_controls==1    
[var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D);
end


%% STEP 6: REPORTING
s=['Results: Leave One Out Largest Connected Set'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
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
s=['Covariance of Firm and Person Effects: ' num2str(sigma_psi_alpha)];
disp(s);
s=['Variance of Person Effects: ' num2str(sigma2_alpha)];
disp(s);
s=['Correlation of Firm Effect and Person Effects: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
disp(s);
s=['Total Explained Variation: ' num2str(explained_var_leave_out)];
disp(s);
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
s=['Leave One Out SE: ' num2str(sqrt(V_sigma2_psi))];
disp(s)
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
s=['Leave One Out SE: ' num2str(sqrt(V_sigma_alpha_psi))];
disp(s)
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
s=['Leave One Out SE: ' num2str(sqrt(V_sigma2_alpha))];
disp(s)
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
%% Focus on Diagnostics
%Note: The eigenvalue ratios (and sum of squared eigenvalues) 
%      reported in Table 3 of KSS have been calculated exactly whereas here
%      we report the results obtained via simulations. 
%      The differences are neglible as one can see from below.
if eigen_diagno==1
for pp=1:Ndiagno
if pp == 1
title='Diagnostics on Variance of Firm Effects';
end    
if pp == 2
title='Diagnostics on Variance of Person Effects';
end
if pp == 3
title='Diagnostics on CoVariance of Person, Firm Effects';
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
for pp=1:Ndiagno
if pp == 1
title='Inference on Variance of Firm Effects';
end    
if pp == 2
title='Inference on Variance of Person Effects';
end
if pp == 3
title='Inference on CoVariance of Person, Firm Effects';
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
%Save File.
s=['mat/COMPLETE_' filename];
save(s)
end

