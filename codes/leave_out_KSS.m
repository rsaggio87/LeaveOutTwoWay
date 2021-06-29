function [sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_KSS(y,id,firmid,controls,leave_out_level,type_algorithm,simulations_JLA,lincom_do,Z_lincom,labels_lincom,filename)
%% Author: Raffaele Saggio
%Email: rsaggio@mail.ubc.ca

%% DESCRIPTION
%This function computes the bias-corrected components in a two-way model as
%described in Kline, Saggio and Soelvsten (2020, ECTA -- KSS Henceforth).

%This function can be applied to any two-way fixed effects model
%(student-teachers, patient-doctors, etc). We use AKM jargon (workers, firms) 
%when describing the code.

%% WARNING!
%This function only works properly if the user has properly sorted the data
%by id-year (e.g. xtset id year in Stata).

%% MANDATORY INPUTS
%y: outcome. Dimensions: N* x 1; N*= # of person-year observations.
%--
%id: worker indicators. Dimensions: N* x 1
%--
%firmid: firm indicators. Dimensions: N* x 1

%% NON-MANDATORY INPUTS    
%leave_out_level: string variable that takes two values:

    %'obs': perform leave-out by leaving a person-year observation out.

    %'matches': perform leave-out by leaving an entire person-firm match out.
    
    %Default: `matches'.

%controls:
    %Matrix of controls with dimensions: N* x P. This matrix of controls must
    %be appropriately defined by the user ex-ante. For instance, if the user 
    %wants to include time effects then the user should include in the matrix 
    %'controls' the set of dummy variables associated to a particular year
    %effect, making sure to avoid potential collinearity issues. 
    
    %These controls are going to be partialled out as follows. Let the
    %model be
    
    %y=D*alpha + F*psi + X*b + e

    %The code is going to estimate the above model, get an estimate of b
    %compute then ynew=y-Xb and run a variance decomposition on 
    
    %ynew=D*alpha + F*psi
    
%type_algorithm: This takes two values: "exact" or "JLA".                    

    %   "exact": performs exact computation of (Bii,Pii). 

    %   "JLA": perform random projection methods to approximate (Bii,Pii) as 
    %   detailed in the Computational Appendix of KSS. 

    %In larger datasets, the user should always set type_algorithm='JLA'.
    
    %Default: JLA if # of observations in the original data is >10,000
    
%simulations_JLA: a natural number.    

    %This governs the # of simulations in the JLA algorithm to approximate 
    %(Bii,Pii).
    
    %Default: 200.
    
%lincom_do: binary. 

    %If =1, the code regresses the firm effects on the set of covariates
    %specied by Z_lincom and report the correct t-statistic
    %associated with this regression.
    
    %Default: 0
    
 %Z_lincom: matrix of regressors with dimension N* x r.
 
    %Matrix of observables to be used when projecting the firm effects into
    %observables.
 
 %labels_lincom: string
 
    %vector of dimension rx1 that provides a label for each of 
    %the columns in Z_lincom.
    
%filename: string. 
    %Where saved results should be stored and named. Use name like 
    %"leave_out_results" and not "leave_out_results.csv" %Default is 'leave_out_estimates';

%% OUTPUTS
    
%sigma2_psi:      Variance of the firm effects.

%sigma_psi_alpha: Covariance of the firm, person effects.

%sigma2_alpha:    Variance of the person effects. 
%                 When leaving a match-out, this parameter is only
%                 identified among movers, i.e. individuals that moved b/w
%                 firms across periods. The log file is going to report bounds
%                 for this parameter + provide its estimate when leaving a
%                 single observation out in the KSS procedure.
%

%The function  is also going to save on disk one .csv file and one .mat file. 
%The .csv contains information on the leave-out connected set. 
%First column reports the outcome variable, second and third columns the 
%worker and the firm identifier. 
%The fourth column reports the stastistical leverages. 
%If the code is reporting a leave-out correction at the match-level, 
%the .csv will be collapsed at the match level. 

%% READ
no_controls=0;
no_algo=0;
no_scale=0;
no_labels=0;
if nargin < 3
    error('More arguments needed');
end

if nargin == 3
    no_controls=1;
    controls=ones(size(y,1),1);
    leave_out_level='matches';
    no_algo=1;
    no_scale=1;
    lincom_do=0;
    no_labels=1;
    filename='leave_out_estimates';
end

if nargin == 4   
    leave_out_level='matches';
    no_algo=1;
    no_scale=1; 
    lincom_do=0;
    no_labels=1;
    filename='leave_out_estimates';
end

if nargin == 5
    no_algo=1;
    no_scale=1;   
    filename='leave_out_estimates';
    lincom_do=0;
    no_labels=1;
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
end

if nargin == 6    
    no_scale=1;   
    filename='leave_out_estimates'; 
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    lincom_do=0;
    no_labels=1;
end

if nargin == 7      
    filename='leave_out_estimates';
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    if size(simulations_JLA,2)==0
    no_scale=1;
    end
    lincom_do=0;
    no_labels=1;
end

if nargin == 8      
    filename='leave_out_estimates';
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    if size(simulations_JLA,2)==0
    no_scale=1;
    end
    if size(lincom_do,2)==0
    lincom_do=0;
    end
    if lincom_do==1 
    disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates')    
    lincom_do=0;
    end
    no_labels=1;
end

if nargin == 9      
    filename='leave_out_estimates';
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    if size(simulations_JLA,2)==0
    no_scale=1;
    end
    if size(lincom_do,2)==0
    lincom_do=0;
    end
    if lincom_do==1 && size(Z_lincom,2)==0 
    disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates,no lincom results will be shown')    
    lincom_do=0;
    end
    no_labels=1;
    
end

if nargin == 10      
    filename='leave_out_estimates';
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    if size(simulations_JLA,2)==0
    no_scale=1;
    end
    if size(lincom_do,2)==0
    lincom_do=0;
    end
    if lincom_do==1 && size(Z_lincom,2)==0 
    disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates,no lincom results will be shown')    
    lincom_do=0;
    end
    if size(labels_lincom,2)==0
    no_labels=1;
    end
    
end

if nargin == 11      
    if size(leave_out_level,2)==0
    leave_out_level='matches';
    end
    if size(type_algorithm,2)==0
    type_algorithm='JLA';
    end
    if size(simulations_JLA,2)==0
    no_scale=1;
    end
    if size(lincom_do,2)==0
    lincom_do=0;
    end
    if lincom_do==1 && size(Z_lincom,2)==0 
    disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates, no lincom results will be shown')    
    lincom_do=0;
    end
    if size(labels_lincom,2)==0
    no_labels=1;
    end
    
end


%Read empty controls
if size(controls,2)==0
    no_controls=1;
    controls=ones(size(y,1),1);
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



%Read # of FE
if no_scale == 1
    simulations_JLA=200; %default rule
end

if no_algo == 1
    
    if size(y,1)>10000
        type_algorithm='JLA';
    end
    
    if size(y,1)<=10000
        type_algorithm='exact';     
    end

end

if 1 == 1
%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp('Running KSS Correction with the following options')
if strcmp(leave_out_level,'matches') 
s=['Leave Out Strategy: Leave match out'];
disp(s);
end
if strcmp(leave_out_level,'obs') 
s=['Leave Out Strategy: Leave person-year observation out'];
disp(s);
end
if strcmp(type_algorithm,'exact')
s=['Algorithm for Computation of Statistical Leverages: Exact'];
disp(s);
end
if strcmp(type_algorithm,'JLA')
s=['Algorithm for Computation of Statistical Leverages: JLA with ' num2str(simulations_JLA) ' simulations.'];
disp(s);
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
end


%% STEP 1: FIND CONNECTED SET
%As first step in our analysis, we run estimation of a standard AKM model
%on the largest connected set

%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Find connected_set
if lincom_do==1
K=size(controls,2);
controls=[controls Z_lincom];
end
[y,id,firmid,id_old,firmid_old,controls] = connected_set(y,id,firmid,lagfirmid,controls);


%% STEP 2: LEAVE ONE OUT CONNECTED SET
%Here we compute the leave out connected set as defined in the Computational
%Appendix of KSS. 
%The input data is represented by the largest connected set. After applying
%the function 'pruning_unbal_v3', the output data will be a connected set
%such that the associated bipartite graph between workers and firms remains
%connected after removing any particular worker from the graph.


%Leave One Out Largest Connected Set
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp('SECTION 1')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
[y,firmid,id,id_old,firmid_old,controls] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);


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
var_den=var(y);
clear id_mover
%Summarize
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
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
%% STEP 3: Residualizing and Collapsing
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp('SECTION 2')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
%If the model includes controls, we're going to estimate the coefficients
%of those controls and partial them out as follows
if lincom_do==1
   Z=controls(:,K+1:end);
   controls=controls(:,1:K); 
end

%Residualize
if no_controls == 0
   S=speye(J-1);
   S=[S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
   X=[D,F*S,controls];
   xx=X'*X;
   xy=X'*y;
   Lchol=lchol_iter(xx);
   if size(Lchol,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
            b           = pcg(xx,xy,1e-10,1000,Lchol,Lchol');
   else 
            b           = pcg(xx,xy,1e-10,1000);
   end
   y=y-X(:,N+J:end)*b(N+J:end); %variance decomposition will be based on this residualized outcome.
end


%Collapsing
%If the user wants to run the KSS correction leaving a match out, we're
%going to collapse the data down to match-means and weight this collapsed
%data by the length of a given match. 

peso=ones(size(y,1),1);
y_py=y;
if strcmp(leave_out_level,'matches') 
   [~,~,match_id] 		= unique([id firmid],'rows','stable');
   peso					= accumarray(match_id,1);
   id					= accumarray(match_id,id,[],@(x)mean(x));
   firmid			    = accumarray(match_id,firmid,[],@(x)mean(x));
   id_old               = accumarray(match_id,id_old,[],@(x)mean(x));
   firmid_old           = accumarray(match_id,firmid_old,[],@(x)mean(x));
   y					= accumarray(match_id,y,[],@(x)mean(x)); 
   if lincom_do==1
   Z					= accumarray(match_id,Z,[],@(x)mean(x));
   end
end

%% STEP 4: COMPUTATION OF (Pii,Bii)
%This is the computationally expensive part of the code where we compute
%the terms (Pii,Bii) as defined in KSS.

%Build Design
    NT=size(y,1);
    D=sparse(1:NT,id',1);
    F=sparse(1:NT,firmid',1);
    X=[D,-F]; %shaped in a pure Laplacian format.
    N=size(D,2);
    J=size(F,2);
    
%Weighting Matrices
	X_fe=[sparse(NT,N) X(:,N+1:end)];
    X_fe=repelem(X_fe,peso,1); %weight by lenght of the spell
    X_fe=[sparse(NT,N) X(:,N+1:end)];
    X_pe=[F*(inv(F'*F)*F'*X(:,1:N)) sparse(NT,J)];
    X_pe=repelem(X_pe,peso,1); %weight by lenght of the spell
    PESO_MAT=sparse(1:NT,(1:NT)',peso.^0.5,NT,NT);
    y_untransformed=y;
    X=PESO_MAT*X;% TO ACCOUNT FOR WEIGHTING (FGLS)
    y=PESO_MAT*y;%FGLS transformation
    xx=X'*X;
    disp('Calculating the statistical leverages...')
    [results,Lchol] = evalc('cmg_sdd(xx)'); %preconditioner for Laplacian matrices.   

tic    
if n_of_parameters==1
    [Pii, Mii, correction_JLA, Bii_fe]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
end

if n_of_parameters==2
    [Pii, Mii, correction_JLA, Bii_fe, Bii_cov]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
end

if n_of_parameters==3
    [Pii, Mii, correction_JLA, Bii_fe, Bii_cov, Bii_pe]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
end

disp('Done!')
toc


%% STEP 5: ESTIMATION OF VARIANCE COMPONENTS
%We use the statistical leverages, Pii, and the Bii associated with a given
%variance component to bias correct these quantities using the KSS approach

%X_fe                = [sparse(NT,N) F];
%X_fe                = repelem(X_fe,peso,1); %weight by lenght of the spell
%X_pe                = [D sparse(NT,J)];
%X_pe                = repelem(X_pe,peso,1); %weight by lenght of the spell
S                   = speye(J-1);
S                   = [S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
X                   = [D,F*S]; %back to grounded Laplacian. 
X                   = PESO_MAT*X;% TO ACCOUNT FOR WEIGHTING (FGLS)
xx                  = X'*X;
xy                  = X'*y;
Lchol               = lchol_iter(xx);
if size(Lchol,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
        b           = pcg(xx,xy,1e-10,1000,Lchol,Lchol');
else 
        b           = pcg(xx,xy,1e-10,1000);
end
eta                 = y-X*b;
eta_h				= eta./Mii; %Leave one out residual
sigma_i				= (y-mean(y)).*eta_h; %KSS estimate of individual variance.
sigma_i				= sigma_i.*correction_JLA; %to adjust for non-linear bias induced by JLA.

X_fe                = X_fe(:,1:end-1);
X_pe                = X_pe(:,1:end-1);
corr                = sigma_i.*Bii_pe;

out=[full(1-Mii),eta_h,full(Bii_pe),full(Bii_fe),corr];
s=['diagno.csv'];
dlmwrite(s, full(out), 'delimiter', '\t', 'precision', 16);


if strcmp(leave_out_level,'matches')
    T               = accumarray(id,1);
    stayers         = T==1;
    stayers         = stayers(id);
    sigma_stayers   = sigma_for_stayers(y_py,id,firmid,peso,b);
    sigma_i(stayers)= sigma_stayers(stayers); 
    
end

if n_of_parameters==1
    [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
end

if n_of_parameters==2
    [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
    [sigma_alpha_psi_AKM, sigma_psi_alpha]  = kss_quadratic_form(sigma_i,X_fe,X_pe,b,Bii_cov);
end

if  n_of_parameters==3
    [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
    [sigma_alpha_psi_AKM, sigma_psi_alpha]  = kss_quadratic_form(sigma_i,X_fe,X_pe,b,Bii_cov);
    [sigma_2_alpha_AKM, sigma2_alpha]       = kss_quadratic_form(sigma_i,X_pe,X_pe,b,Bii_pe);
end

    pee=b(1:N);
    pee=X(:,1:N)*pee;
    pee=accumarray(firmid,pee,[],@(x)mean(x));
    var(pee)

%% STEP 6: PRINTING RESULTS
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
disp('SECTION 3')
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
s=['PLUG-IN ESTIMATES (BIASED)'];
disp(s)
s=['Variance of Firm Effects: ' num2str(sigma_2_psi_AKM)];
disp(s)
if n_of_parameters>=2
s=['Covariance of Firm, Person Effects: ' num2str(sigma_alpha_psi_AKM)];
disp(s);
end
if n_of_parameters==3
s=['Variance of Person Effects: ' num2str(sigma_2_alpha_AKM)];
disp(s);
s=['Correlation of Firm, Person Effects: ' num2str(sigma_alpha_psi_AKM/(sqrt(sigma_2_psi_AKM)*sqrt(sigma_2_alpha_AKM)))];
disp(s);
s=['Fraction of Variance explained by Worker and Firm Effects: ' num2str((sigma_2_psi_AKM+2*sigma_alpha_psi_AKM+sigma_2_alpha_AKM)/var_den)];
disp(s);
end
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
s=['BIAS CORRECTED ESTIMATES'];
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
s=['Correlation of Firm, Person Effects: ' num2str(sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)))];
disp(s);
s=['Fraction of Variance explained by Worker and Firm Effects: ' num2str((sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var_den)];
disp(s);
end
%% STEP 7: SAVING OUTPUT
%Export csv with data from leave out connected set
out=[y_untransformed,id_old,firmid_old,full(Pii)];
s=[filename '.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
%% STEP 8: LINCOM
if lincom_do == 1 
    Z                   = repelem(Z,peso,1);
    Transform           = X_fe; %regress the firm effects (+ regression is person-year weighted)
    disp('Regressing the firm effects on observables...')
if no_labels == 0
    [~]                 = lincom_KSS(y,X,Z,Transform,sigma_i,labels_lincom);
end
if no_labels == 1
    [~]                 = lincom_KSS(y,X,Z,Transform,sigma_i);
end
end
end