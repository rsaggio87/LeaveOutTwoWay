function tabella = leave_out_KSS_coll(y,id,firmid,Tspell,leave_out_level,type_algorithm,simulations_JLA,filename,group_indicator,group_index)

if 1 == 1
%Listing options
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp('Running KSS Correction on data collapsed at the match level')
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
controls=[group_indicator Tspell];
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


%% Dropping stayers from the analysis
disp('Dropping stayers....')
T=accumarray(id,1);
T=T(id);
sel=T>1;
y=y(sel,:);
firmid=firmid(sel,:);
id=id(sel,:);
id_old=id_old(sel,:);
firmid_old=firmid_old(sel,:);
controls=controls(sel,:);
group_indicator=controls(:,1);
peso=controls(:,2);
%Resetting ids one last time.
[~,~,n]=unique(firmid);
firmid=n;
[~,~,n]=unique(id);
id=n; 

%Important Auxiliaries
conteggio_matches           = sum(group_indicator==group_index);
y_expand                    = repelem(y,peso,1);
group_indicator_coll        = group_indicator;
group_indicator             = repelem(group_indicator,peso,1);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Summary-Stats-->group:' num2str(group_index) ' (KSS leave-out connected set; all stats are person-year weighted)'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
sel = group_indicator == group_index;
conteggio_py                = sum(sel);

%% py weighted stats
s=['mean outcome: ' num2str(mean(y_expand(sel)))];
disp(s)
s=['variance of match means: (might want to adjust this so that it includes also within-job variability) ' num2str(var(y_expand(sel)))];
disp(s)
mean_outcome=mean(y_expand(sel));
var_den=var(y_expand(sel));
y_expand=[];
%% count workers and firms in group g
sel = group_indicator_coll == group_index;
id_sel = id(sel);
[~,~,id_sel]=unique(id_sel);
s=['# of Workers: ' num2str(max(id_sel))];
disp(s);
id_sel = firmid(sel);
[~,~,firmid_sel]=unique(id_sel);
s=['# of Firms: ' num2str(max(firmid_sel))];
disp(s);
s=['# of Person Year Observations: ' num2str(conteggio_py)];
disp(s);
s=['# of Matches: ' num2str(conteggio_matches)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
N_of_workers    = max(id_sel);
N_of_firms      = max(firmid_sel);
N_of_OBS        = sum(sel);
var_outcome     = var_den;


if strcmp(leave_out_level,'matches') 
   disp('error: you cannot run leave-out matches on an already collapsed data')
end

%% STEP 3: COMPUTATION OF (Pii,Bii)
%This is the computationally expensive part of the code where we compute
%the terms (Pii,Bii) as defined in KSS.

%Build Design Matrices
    NT=size(y,1);
    D=sparse(1:NT,id',1);
    F=sparse(1:NT,firmid',1);
    J = max(firmid);
    S=speye(J-1);
    S=[S;sparse(-zeros(1,J-1))];
    X=[D,F*S]; %shaped in a pure Laplacian format.
    N=size(D,2);
    J=size(F,2);
    
%Weighting Matrices for Variance Components
	X_fe=[sparse(NT,N) X(:,N+1:end)];
    X_fe=repelem(X_fe,peso,1); %weight by lenght of the spell
    sel = group_indicator == group_index;
    X_fe=X_fe(sel,:);
    X_pe=[X(:,1:N) sparse(NT,J-1)];
    X_pe=repelem(X_pe,peso,1); %weight by lenght of the spell
    X_pe=X_pe(sel,:);
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
    disp('Done!')
    toc


%% STEP 4: ESTIMATION OF VARIANCE COMPONENTS
%We use the statistical leverages, Pii, and the Bii associated with a given
%variance component to bias correct these quantities using the KSS approach
    X_fe                = [sparse(NT,N) F];
    X_fe                = repelem(X_fe,peso,1); %weight by lenght of the spell
    X_fe                = X_fe(sel,:);
    X_pe                = [D sparse(NT,J)];
    X_pe                = repelem(X_pe,peso,1); %weight by lenght of the spell
    X_pe                = X_pe(sel,:);
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

%% STEP 6: PRINTING RESULTS
    s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
    disp(s);
    disp(s);
    s=['Variance Decomposition for group: ' num2str(group_index)];
    disp(s);
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
    s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
    disp(s);
    disp(s);
    end

%% STEP 7: SAVE OUTPUT
%Export csv with data from leave out connected set with person, firm
%effects
cs      = [id_old firmid_old fe pe id firmid];
cs      = full(cs);
s=['results/' filename '.csv'];
dlmwrite(s, cs, 'delimiter', '\t', 'precision', 16);

%Export variance decomposition within each subgroup
tabella = [N_of_workers; N_of_firms; N_of_OBS; -9;-9;mean_outcome; var_outcome;-9;-9; sqrt(sigma2_psi); sqrt(sigma2_alpha);  sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha)); -9; -9;-9; sigma2_psi/var_den; sigma2_alpha/var_den; (2*sigma_psi_alpha)/var_den];
tabella=full(tabella);
s=['results/variance_decomp___GROUP' num2str(group_index) '___' filename '.csv'];
dlmwrite(s, tabella, 'delimiter', '\t', 'precision', 16);

