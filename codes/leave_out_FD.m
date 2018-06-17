function [sigma2_psi,V] = leave_out_FD(y,id,firmid,leave_out_level,controls,type_algorithm,eigen_diagno,eigen_fast,do_montecarlo,filename)
%% Author: Raffaele Saggio
%Email: raffaele.saggio@berkeley.edu

%% Version:
% 1.0: Wrote documentation. 06.15.2018.


%% DESCRIPTION
%This function computes the leave out estimates of the variance
%of firm effects in two-way fixed effects model as described in KSS.

%This functions permits computation of the leave-out variance of firm 
%effects also in very large datasets by partialling out the worker effects
%and working with the model in first differences (or in stacked First
%Differences when T>2, see below).

%When working with large dataset it is highly recommended to set the option
%type_algorithm='JL', in order to apply the randomized algorithm routine 
%described in Appendix B of KSS. 
%Also, the user should set the option eigen_fast=1 in order to speed up 
%calculations of the eigenvalues and eigevectors.

%When the input data set has max(T_i)=2, where T_i is the total number of
%person year observations in which we observe a worker, the function works
%with a simple model in First Differences (FD). When max(T_i)>2, we work
%with a model that combines all first differences for a given individual and
%notice that weighted least squares estimates in this Pooled First Differences
%(PFD) model are numerically equivalent to standard within group 
%(Fixed Effects) estimates, provided that each difference for an individual
%is weighted by T_i.

%The code also includes a small montecarlo exercise at the end as described
%in the empirical section of KSS.

%If controls are specifiied by the user the function will work as follows:
%After finding the largest connected set the code estimates a standard AKM
%model with controls. Then it uses the estimated effects on these controls to 
%partial out the effect of these variables on the outcomes. The leave-out
%model will then work with this residualized outcome to speed up
%computation.

%Note: future releases of the function 'leave_out_complete' will also allow
%fast computation of both the variance of person effects and the 
%covariance of person, firm effects, see https://sites.google.com/site/raffaelesaggio/


%% DESCRIPTION OF THE INPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.

%id: worker indicators. Dimensions: N* x 1

%firmid: firm indicators. Dimensions: N* x 1

%leave_out_level: string variable that takes three values:

%'obs': perform leave-out by leaving a person-year observation out (default)

%'workers': perform leave-out by leaving an entire worker's history out.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                    %---NON-MANDATORY INPUTS
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%type_algorithm: 'string'

%This declares which type of algorithm to perform in order to
%compute (Bii,Pii). 

%If type_algorithm='exact', then (Bii,Pii)
%are computed using an exact method which can take quite a bit of time in 
%large datasets. 

%If type_algorithm='JL', then the code uses the Spielman and Srivastava (2008) 
%algorithm which make use of the Johnson-Lindenstrauss lemma. 
%See Appendix B

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%controls:
%Matrix of controls with dimensions: N* x K. This matrix of controls must
%be appropriately defined by the user ex-ante. For instance, if the user 
%wants to include time effects then the user should include in the matrix 
%'controls' the set of dummy variables associated to a particular year
%effect. If 'controls' is empty, then no controls will be used for
%estimation.


%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                  
%eigen_diagno: Binary. 

% If 1, the code outputs the lindeberg condition and 
% eigenvalue ratio of theorem 1. The code will also output the 
% weak-id confidence intervals using the AM method described in the paper 
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%eigen_fast: Binary. 
%If 1 (and eigen_diagno=1), the code uses a simulation method to estimate 
%the sum of squared eigenvalues of the design matrix.
%It is recommended to set eigen_fast=1 when working with large datasets.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%do_montecarlo: Binary. 
%If 1, the code performs a MC experiment.
% Default is 0.

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%filename: string. 
%Used to name saved outputs.
%Default is 'leave_out_FD_estimates';
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-


%% DESCRIPTION OF THE OUTPUTS

%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
%sigma2_psi: leave-out variance of firm effects. 
%V:sampling variance of the leave-out variance of firm effects

%Check Log File for additional results reported by the code. 
%-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

%% EXTERNAL CALLS
MakeCMG; %preconditioner for Laplacians matrix: http://www.cs.cmu.edu/~jkoutis/cmg.html
path(path,'~/matlab_bgl/'); %path to the matlabBGL files. note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922

no_controls=0;
%% READ
if nargin < 4
error('More arguments needed');
end

if nargin == 4
    no_controls=1;
    controls=ones(size(y,1),1);
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
end

if nargin == 5
    type_algorithm='JL';
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
end

if nargin == 6
    eigen_diagno=0;
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 7
    eigen_fast=0;
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 8
    do_montecarlo=0;
    filename='leave_out_FD_estimates';
end

if nargin == 9
    filename='leave_out_FD_estimates';
end

if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
end



%% Identify the type of panel (T=2 or T>2).
[~,~,id_norm]=unique(id);
T=accumarray(id_norm,1);
clear id_norm
maxT=max(T);

if maxT==2 % T=2.
    type_estimator='FD';
end

if maxT>2 % General case. Will set model in PFD
    type_estimator='FE';
end


if do_montecarlo==1 && strcmp(type_estimator,'FE')  
    error('Montecarlo can be computed only for T=2 case.');
end



%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
disp('Listing options')
no_controls
leave_out_level
type_algorithm
type_estimator
eigen_diagno
eigen_fast
do_montecarlo
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
L=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
toc
b=pcg(xx,xy,1e-10,1000,L,L');
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
%Do some cleaning of matrices in memory
clear xx xy L xb pe fe ahat ghat F D S Lchol


%%Deresidualized outcome variable (allows the leave one out methodology to run faster)
if no_controls==0
y=y-X(:,N+J:end)*b(N+J:end);
end
clear X b


%% STEP 2: LEAVE ONE OUT CONNECTED SET
%Here we compute the leave out connected set as defined in Appendix B. 
%The input data is represented by the largest connected set. After applying
%the function 'pruning_unbal_v3', the output data will be a connected set
%such that the associated bipartite graph between workers and firms remains
%connected after removing any particular worker from the graph.

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


%Very Important Auxiliaries
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
clear id_mover

s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
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

%% STEP 3: RESHAPING DATA
%If max(T_i)=2, we simply set the model in FD.
%If max(T_i)>2, we set the model in PFD and use the appropriate weighting.

NT=size(y,1);
F=sparse(1:NT,firmid',1);
J=size(F,2);
count=ones(length(id),1);
gcs = cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));


%Standard FD transformation
if strcmp(type_estimator,'FD')
    %outcome variable    
    ylag=[NaN; y(1:end-1)];
    ylag(gcs==1)=NaN; %%first obs for each worker
    ydelta=y-ylag;
    ydelta(isnan(ydelta(:,1)),:)=[]; %remove first observation for worker;

    %matrix of assignements
    Flag=[NaN(1,J); F(1:end-1,:)];
    Flag(gcs==1,1)=NaN; %just put it for first firm, we'll remove rows where there is at least one NaN.
    Fdelta=F-Flag;
    sel=~any(isnan(Fdelta),2);
    Fdelta=Fdelta(sel,:);

    %ids
    firmid_delta=[NaN; firmid(1:end-1)];
    firmid_delta(gcs==1)=NaN;
    firmid_delta=firmid_delta(sel);
    firmid_delta_f=firmid;
    firmid_delta_f=firmid_delta_f(sel);
    id_movers=id(sel);
    clear Flag ylag
end

%PD transformation
if strcmp(type_estimator,'FE')
    tic
    [ydelta, Fdelta, Tweight,id_movers,firmid_delta,firmid_delta_f,gcs_delta]= stacked_Fdelta(y,id,firmid,gcs);
    disp('Time to set up the model in Pooled Differences')
    toc
end

%Design matrix (Laplacian) and preconditioner
L=(Fdelta'*Fdelta); %Laplacian matrix
pfun_ = cmg_sdd(L); %preconditioner for Laplacian matrices.

%% STEP 4: DIAGNOSTICS FOR ESTIMATING THE VARIANCE OF FIRM EFFECTS
%This part is only computed provided that the option 'eigen_diagno' is
%turned on. In this part, we calculate the squared eigenvalue ratio and
%Lindeberg conditions that constitute the key conditions to verify the
%validity of Theorem 1. Lindebergc condition assumes q=1.

if eigen_diagno==1
    EIG_NORM=zeros(3,1);
    [Q, lambda_eig] = eigs(F'*F,L,4);
    lambda_eig=diag(lambda_eig);
    lambda_eig=lambda_eig(2:end); %Laplacian has always first eigenvalue equal to 0.
    lambda_1=lambda_eig(1);
    q=Q(:,1);
    x1bar=Fdelta*q;
    norm=(sum(x1bar.^2))^(0.5);
    x1bar=x1bar/norm;
    disp('checking, must report 1')
    sum(x1bar.^2)
    clear Q
    
    if eigen_fast==0 %To calculate sum of squared eigenvalues (here using exact method)
        Degree=(F'*F);
        Dsqrt=sqrt(Degree);
        Dsqrt_inv=Dsqrt^(-1);
        normL=Dsqrt_inv*L*Dsqrt_inv;
        EIG = eig(full(normL));
        EIG = sort(EIG);
        EIG=EIG(2:end); %Laplacian has always first eigenvalue equal to 0.
        EIG=(1./EIG); %eigenvalues for the inverse.
        EIG=EIG.^2;
        SUM_EIG=sum(EIG);
        clear normL Degree Dsqrt_inv Dsqrt
            for pp=1:3
                EIG_NORM(pp,1)=EIG(pp)/SUM_EIG;
            end
    end
    
    if eigen_fast==1 %To calculate sum of squared eigenvalues (via simulations)
        SUM_EIG=trace_Atilde_sqr_FD(Fdelta,F,L,pfun_); 
            for pp=1:3
                EIG_NORM(pp,1)=(lambda_eig(pp,1)^2)/SUM_EIG;
            end    
    end
    disp('Time to find Eigenvalues and Eigenvectors to compute Lindeberg')
    toc
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%% STEP 4: IDENTIFY THE NON ZERO VALUES IN HAT MATRIX
%The next step is probably the most abtstract, yet the most important,
%passage of the code. Leave one out computation requires computation of the
%(block) diagonal of the "hat" matrix, P. This is potentially a
%gigantic matrix. Even if it is block diagonal, some computers might not be
%able to store all its elements. However, thanks to the particular
%structure associated with the Laplacian we know exactly which entries 
%need to be filled in this matrix. 

%If the user perfoms leave-out dropping just a person-year observation,
%then the associated (Bii,Pii) will be a non-zero scalar if the worker
%moved across firms in the corresponding year-year combination. 

%If the user perfoms leave-out dropping the entire history of a worker,
%then the associated (Bii,Pii) will each be a matrix as we have to keep
%track of cross-products across year-year combinations where the worker has
%moved across firms.

tic
if strcmp(type_estimator,'FD')
    index=[(1:size(firmid_delta,1))']; %index of observations in the FD world
    sel=firmid_delta~=firmid_delta_f;
    elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since here T=2, no need for cross-products when doing leave-out
    rows=index(sel);
    column=index(sel);
    [elist, ~, index_unique]=unique(elist,'rows');
end

if strcmp(type_estimator,'FE') && strcmp(leave_out_level,'obs') 
    index=[(1:size(firmid_delta,1))']; %index of observations in the PFD world
    sel=firmid_delta~=firmid_delta_f;
    elist=[firmid_delta(sel) firmid_delta_f(sel) firmid_delta(sel) firmid_delta_f(sel)]; %keeping track of firm in t, firm in t' for movers. Since we want to leave one obs at the time, no need for cross-products when doing leave-out
    rows=index(sel);
    column=index(sel);
    Tweight=Tweight(sel);
    [elist, ~, index_unique]=unique(elist,'rows');
end

if strcmp(type_estimator,'FE') && strcmp(leave_out_level,'workers')
    [elist,index_unique,rows,column,Tweight]= input_lambda_P_FE_fast(firmid_delta,firmid_delta_f,id_movers,Tweight); %%keeping track of firm in t, firm in t' for movers across periods for a given worker. Now we need also to keep track of cross-products because we leave want to run leave-out leaving the entire worker history of a given individual
end

disp('Time to set up the indexes of the Hat matrix')
toc
disp('Completed Step 4')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%% STEP 5: COMPUTE (Bii,Pii)
%There are two ways to proceed. 

%First way, is to compute exact estimates,
%by parallelizing computation of S_xx^(-1)x_i' across cores.

%Second way uses the SS (2011) algorithm described in Appendix B. Notice
%that everything applies of that Appendix because we still have S_xx has a
%Laplacian matrix in this context.

%We make use of the cmg solver for both steps. The user can fix a seed for 
%replication purposes. 


%Inputs
tol=1e-5; %tol for pcg
epsilon=0.01; %rules the tol for random projection (essentially how many simulations to take).

%Calculate (Bii,Pii)
[Pii, Bii] = eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_);
Pii=Pii(index_unique);
Bii=Bii(index_unique);

%need reweight differences if running FE via PD.
if strcmp(type_estimator,'FE')
    Pii=Pii./Tweight; 
    Bii=Bii./Tweight;
end

%Lambda P
Lambda_P=sparse(rows,column,Pii,size(Fdelta,1),size(Fdelta,1));
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.
%Lambda B
Lambda_B=sparse(rows,column,Bii,size(Fdelta,1),size(Fdelta,1));
Lambda_B=Lambda_B+triu(Lambda_B,1)'; %make it symmetric.
disp('Completed Step 5')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
%save(s)
%% STEP 6: ESTIMATION
%Step 6A: AKM (Plug-in) estimates of the variance of firm effects.
%Step 6B: Verify that AKM and PDF give back same answer.
%Step 6C: Leave out with associated confidence interval (q=0 and q=1).

%% STEP 6A: AKM ESTIMATION on LEAVE ONE OUT LARGEST CONNECTED SET
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix
D=sparse(1:NT,id',1);
X=[D,F*S];
xx=X'*X;
xy=X'*y;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
b=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
N=size(D,2);
J=size(F,2);
ahat=b(1:N);
ghat=b(N+1:N+J-1);
pe=D*ahat;
fe=F*S*ghat;
COV=cov(fe,pe);
sigma_2_psi_AKM=COV(1,1);
clear X b xx Lchol ahat ghat S 
%% STEP 6B: VERIFY model in PFD gives back same estimates as Fixed Effects.
xy=Fdelta'*ydelta;
b=pcg(L,xy,1e-10,1000,pfun_);
psi_hat_norm=b;
disp('checking - must report 0')
abs(var(fe)-var(F*b)) %you can also check that each single firm effect (when consistently normalized across methods) are numerically equivalent.
r=ydelta-Fdelta*b;
clear xy
%% STEP 6C: LEAVE OUT ESTIMATES

%Leave Out Residuals
Ndelta=size(Fdelta,1);
I_Lambda_P=(speye(Ndelta,Ndelta)-Lambda_P);
L_P=ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[eta_h, flag]=pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');

%Auxiliary
A_b=F'*(fe-mean(fe));
[W_to_use, my_first_part] = construc_W_FD(ydelta,Fdelta,L,pfun_,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);

tic
if nargout==1
    sigma2_psi=leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,sigma_2_psi_AKM,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part);
end

if eigen_diagno == 0 && nargout==2
    [sigma2_psi, V, sigma_predict]= leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,sigma_2_psi_AKM,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part);
end

if eigen_diagno == 1 && nargout==2
    [sigma2_psi, V, sigma_predict, COV_R1, gamma_sq,Fstatistic,b_1,theta_1]= leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,sigma_2_psi_AKM,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part,x1bar,lambda_1);
end
disp('Time to Compute Leave Out of Firm Effects')
toc

%% STEP 7: REPORTING
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
s=['Max Leverage' num2str(max(diag(Lambda_P)))];
disp(s);
s=['-*-*-*-*-*-*Variance of Firm Effects-*-*-*-*-*-*'];
disp(s)
s=['AKM: ' num2str(sigma_2_psi_AKM)];
disp(s)
s=['Leave One Out: ' num2str(sigma2_psi)];
disp(s)
if nargout==2
s=['Leave One Out SE: ' num2str(sqrt(V))];
disp(s)
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
if eigen_diagno==1
for pp=1:1
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
s=['Weak Lindeberg Condition: ' num2str(max(x1bar.^2))];
disp(s)
s=['Sum of squared eigenvalues: ' num2str(SUM_EIG/NT^2)];
disp(s)
s=['Variance of b_1:  ' num2str(COV_R1(1,1))];
disp(s);
s=['Variance of theta_1: ' num2str(COV_R1(2,2))];
disp(s);
s=['Covariance of (b_1,theta_1): ' num2str(COV_R1(1,2))];
disp(s);
s=['Correlation of (b_1,theta_1): ' num2str((COV_R1(1,2))/(sqrt(COV_R1(2,2))*sqrt(COV_R1(1,1))))];
disp(s);
s=['gamma squared: ' num2str(gamma_sq)];
disp(s);
s=['Fstatistic: ' num2str(Fstatistic)];
disp(s);
end
end
s=['******************************************'];
disp(s);
%% Focus on Inference
if nargout==2
for pp=1:1
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
s=['SE under q=0: ' num2str(sqrt(V))];
disp(s);
s=['CI under q=0: ' num2str(sigma2_psi-1.96*sqrt(V)) '  '  '  ' num2str(sigma2_psi+1.96*sqrt(V))];
disp(s);   
if eigen_diagno==1
[UB,LB,C]=AM_CI(NT,lambda_1,gamma_sq,COV_R1,b_1,theta_1);
s=['CI under q=1: ' num2str(LB) '  '  '  ' num2str(UB)];
disp(s); 
s=['Curvature: ' num2str(C)];
disp(s); 
s=['******************************************'];
disp(s); 
end
end
end

%Save File.
s=['mat/FD_' filename];
save(s)


%% STEP 8: MONTECARLO
%Run a parametric bootstrap simulation setting the truth as the observed
%network structure. Only works for T=2. We only care about variance of firm effects. 

if do_montecarlo==1 
SIMUL=10000; %number of simulations
dof_t=5;
var_t=dof_t/(dof_t-2);
sd_t=sqrt(var_t);

%True parameters, calibrated according to estimates obtained in original data.
naive_var=var(F*psi_hat_norm);
psi_true=psi_hat_norm*sqrt(sigma2_psi/naive_var);
sigma2_true=var(F*psi_true);
   

%Auxiliaries
NT=size(y,1);
trace_B=trace(Lambda_B);
sel=sigma_predict<0;
sigma_aux=abs(sigma_predict(sel));
sigma_predict(sel)=sigma_aux;
dof=size(Fdelta,1)-J+1;
clear sigma_aux

% Fill up the vectors
sigma2_psi_simul=zeros(SIMUL,1);
sigma2_psi_simul_AKM=zeros(SIMUL,1);
sigma2_psi_simul_andrews=zeros(SIMUL,1);
V_exact_simul=zeros(SIMUL,1);
C_simul=zeros(SIMUL,1);
coverage_exact_95=zeros(SIMUL,1);
coverage_fancy_95=zeros(SIMUL,1);
COV_R1_simul=zeros(2,2,SIMUL);
bar_Beta_1_simul=zeros(SIMUL,1);
theta_2_simul=zeros(SIMUL,1);


parfor s=1:SIMUL
%DGP
eta=(sqrt(sigma_predict).*random('t',dof_t,size(Fdelta,1),1))/sd_t;
ydelta=Fdelta*psi_true+eta;

%AKM
xy=Fdelta'*ydelta;
[psi_hat, flag]=pcg(L,xy,1e-10,1000,pfun_);
fe=F*psi_hat;
r=ydelta-Fdelta*psi_hat;
sigma_2_psi_AKM=var(fe);

%Leave-out
[eta_h, flag]=pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');
A_b=F'*(fe-mean(fe));
[W_to_use, my_first_part] = construc_W_FD(ydelta,Fdelta,L,pfun_,A_b,Lambda_B,I_Lambda_P,L_P,eta_h);
[sigma2_psi_s, V_exact_s, ~, COV_R1, gamma_sq,~,b_1,theta_1]= leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,sigma_2_psi_AKM,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part,x1bar,lambda_1);   

%ANDREWS
MSE=sum(r.^2)/dof;
sigma2_andrews=sigma_2_psi_AKM-(1/(NT))*MSE*trace_B;

%Assign results from simulation
sigma2_psi_simul(s)=sigma2_psi_s;
sigma2_psi_simul_AKM(s)=sigma_2_psi_AKM;
sigma2_psi_simul_andrews(s)=sigma2_andrews;
COV_R1_simul(:,:,s)=COV_R1;
bar_Beta_1_simul(s)=b_1;
theta_2_simul(s)=theta_1;

%SEs
V_exact_simul(s)=V_exact_s;
V_exact_simul(s)=sqrt(V_exact_simul(s));

%Inference at 95%
UB=sigma2_psi_simul(s)+1.96*V_exact_simul(s);
LB=sigma2_psi_simul(s)-1.96*V_exact_simul(s);
coverage_exact_95(s)=(sigma2_true>=LB && sigma2_true<= UB);
[UB,LB,C_simul(s)]=AM_CI(NT,lambda_1,gamma_sq,COV_R1,b_1,theta_1);
coverage_fancy_95(s)=(sigma2_true>=LB && sigma2_true<= UB);
end

%% STEP 11: OUTPUT MONTECARLO

%Oracle Estimator
oracle_UB=sigma2_psi_simul+1.96*std(sigma2_psi_simul);
oracle_LB=sigma2_psi_simul-1.96*std(sigma2_psi_simul);
oracle_coverage=(sigma2_true>=oracle_LB).*(sigma2_true<=oracle_UB);

%Oracle Estimator fancy
COV_R1_oracle=cov(bar_Beta_1_simul,theta_2_simul);
gamma_sq_oracle=((lambda_1^2/NT^2)*(COV_R1_oracle(1,1)^2))/(COV_R1_oracle(2,2));
oracle_coverage_fancy=zeros(SIMUL,1);
parfor s=1:SIMUL
    [UB,LB] = AM_CI(NT,lambda_1,gamma_sq_oracle,COV_R1_oracle,bar_Beta_1_simul(s),theta_2_simul(s))
    oracle_coverage_fancy(s)=(sigma2_true>=LB).*(sigma2_true<=UB);
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['True Value of Variance of Firm Effects: ' num2str(sigma2_true)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of Leave one Out Variance Estimator: ' num2str(mean(sigma2_psi_simul))];
disp(s);
s=['Std of Estimated Variance of Firm Effects: ' num2str(std(sigma2_psi_simul))];
disp(s);
s=['Expected Value of AKM: ' num2str(mean(sigma2_psi_simul_AKM))];
disp(s);
s=['Std Deviation  of AKM: ' num2str(std(sigma2_psi_simul_AKM))];
disp(s);
s=['Expected Value of Andrews: ' num2str(mean(sigma2_psi_simul_andrews))];
disp(s);
s=['Std Deviation  of Andrews: ' num2str(std(sigma2_psi_simul_andrews))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of standard error estimator (exact when q=0): ' num2str(mean(V_exact_simul))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Coverage Rate of oracle estimator 95%: ' num2str(mean(oracle_coverage))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (exact) 95%: ' num2str(mean(coverage_exact_95))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (AM Method) 95%: ' num2str(mean(coverage_fancy_95))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (AM Method - oracle) 95%: ' num2str(mean(oracle_coverage_fancy))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['p1 p10 p25 p50 p75 p90 p99 of mathcal{C}%:']
disp(s)
quantile(C_simul,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)


% Oracle Standardized
oracle_norm=(sigma2_psi_simul-sigma2_true)./(std(sigma2_psi_simul));

% Export File for Plots of Simulation Results
s=['mat/oracle_norm' filename '.csv'];
out=oracle_norm;
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 

% Export File for Matlab
s=['/mat/Montecarlo_' filename];
save(s)
end 
end

