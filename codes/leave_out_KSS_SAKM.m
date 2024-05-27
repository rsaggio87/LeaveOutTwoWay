function [cov_person_effects_dual var_person_effects_dual_union var_person_effects_dual_non_union] = leave_out_KSS_SAKM(y,id,firmid,id_orig,controls,leave_out_level,type_algorithm,simulations_JLA,filename,union_status,year)

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
n_of_parameters=3;
%% STEP 1: FIND CONNECTED SET
%As first step in our analysis, we run estimation of a standard AKM model
%on the largest connected set

%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker



%Find connected_set
K=size(controls,2);
controls=[controls year id_orig union_status];
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
[results,y,firmid,id,id_old,firmid_old,controls] = evalc('pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls)');


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
union_status=controls(:,end);
id_orig=controls(:,end-1);
year=controls(:,end-2);
controls=controls(:,1:end-3);

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

%Calculate number of workers who moved their union status
[~,~,n]=unique(id_orig);
id_orig_norm=n; 
switched = accumarray(id_orig_norm,union_status,[],@(x)mean(x));
switch_union = 1.*(switched>0 & switched<1);
switch_union = switch_union(id_orig_norm);

clear id_mover
%If the model includes controls, we're going to estimate the coefficients
%of those controls and partial them out as follows
%Summarize
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Info on the leave one out connected set'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
var_den = var(y);
disp(s)
s=['# of Workers: ' num2str(max(id_orig_norm))];
disp(s);
s=['# of Workers by Union: ' num2str(max(id))];
disp(s);
s=['# of Firms: ' num2str(max(firmid))];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['Share of Switchers (p-y weighted): ' num2str(mean(switch_union))];
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
%Residualize
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
   cs=[id_old year firmid_old id_orig union_status y];

%Matrices for switchers
  [X_union, X_non_union]=mover_union(firmid,id,id_orig,union_status);
    
   
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
   y					= accumarray(match_id,y,[],@(x)mean(x));
end

%% STEP 4: COMPUTATION OF (Pii,Bii)
%This is the computationally expensive part of the code where we compute
%the terms (Pii,Bii) as defined in KSS.

%Build Design
    NT=size(y,1);
    D=sparse(1:NT,id',1);
    F=sparse(1:NT,firmid',1);
    S=speye(J-1);
    S=[S;sparse(-zeros(1,J-1))];
    X=[D,F*S]; %shaped in a pure Laplacian format.
    N=size(D,2);
    J=size(F,2);
    
%Weighting Matrices
	X_fe=[sparse(NT,N) X(:,N+1:end)];
    X_fe=repelem(X_fe,peso,1); %weight by lenght of the spell
    X_pe=[X(:,1:N) sparse(NT,J-1)];
    X_pe=repelem(X_pe,peso,1); %weight by lenght of the spell
    PESO_MAT=sparse(1:NT,(1:NT)',peso.^0.5,NT,NT);
    X=PESO_MAT*X;% TO ACCOUNT FOR WEIGHTING (FGLS)
    y=PESO_MAT*y;%FGLS transformation
    xx=X'*X;
    disp('Calculating the statistical leverages...')
    Lchol               = lchol_iter(xx); %preconditioner for Laplacian matrices.   
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

 %Matrices for dual connected firms
 [~, ~, ~, Bii_union, Bii_cov_union, Bii_non_union]= leverages(X_union,X_non_union,X,xx,Lchol,type_algorithm,simulations_JLA);

disp('Done!')
toc


%% STEP 5: ESTIMATION OF VARIANCE COMPONENTS
%We use the statistical leverages, Pii, and the Bii associated with a given
%variance component to bias correct these quantities using the KSS approach

X_fe                = [sparse(NT,N) F];
X_fe                = repelem(X_fe,peso,1); %weight by lenght of the spell
X_pe                = [D sparse(NT,J)];
X_pe                = repelem(X_pe,peso,1); %weight by lenght of the spell
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
pe                  = b(1:N);
pe                  = D*pe;
fe                  = b(N+1:end);
fe                  = [fe;0];
fe                  = F*fe;

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

[var_person_effects_dual_union_PI, var_person_effects_dual_union]                = kss_quadratic_form(sigma_i,X_union,X_union,b,Bii_union);
[var_person_effects_dual_non_union_PI, var_person_effects_dual_non_union]        = kss_quadratic_form(sigma_i,X_non_union,X_non_union,b,Bii_non_union);
[cov_person_effects_dual_PI, cov_person_effects_dual]                            = kss_quadratic_form(sigma_i,X_union,X_non_union,b,Bii_cov_union);

%% STEP 6: PRINTING RESULTS
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
s=['Variance of Person Effects ---SWITCHERS, WHILE ON UNION: ' num2str(var_person_effects_dual_union_PI)];
disp(s);
s=['Variance of Person Effects ---SWITCHERS, WHILE ON NON-UNION: ' num2str(var_person_effects_dual_non_union_PI)];
disp(s);
s=['SLOPE ---SWITCHERS" ' num2str(cov_person_effects_dual_PI / var_person_effects_dual_union_PI )];
disp(s);
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
R2=(sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var_den;
end
s=['Variance of Person Effects ---SWITCHERS, WHILE ON UNION: ' num2str(var_person_effects_dual_union)];
disp(s);
s=['Variance of Person Effects ---SWITCHERS, WHILE ON NON-UNION: ' num2str(var_person_effects_dual_non_union)];
disp(s);
s=['SLOPE ---SWITCHERS: ' num2str(cov_person_effects_dual/var_person_effects_dual_union)];
disp(s);

%% STEP 7: SAVE OUTPUT
%Export csv with data from leave out connected set with person, firm
%effects and statistical leverages
fe      = repelem(fe,peso,1); %wback to person-year space
pe      = repelem(pe,peso,1); %wback to person-year space
firmid  = repelem(firmid,peso,1); %wback to person-year space
id      = repelem(id,peso,1);
cs      = [cs fe pe id firmid cov_person_effects_dual.*ones(size(fe,1),1) var_person_effects_dual_union.*ones(size(fe,1),1) var_person_effects_dual_non_union.*ones(size(fe,1),1)];
cs      = full(cs);
s=['data/' filename '.csv'];
dlmwrite(s, cs, 'delimiter', '\t', 'precision', 16);
end