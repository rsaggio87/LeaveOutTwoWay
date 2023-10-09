function [coeff,SE,SE_naive]  = KSS_SE(y,D,clusterID,year,controls,labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECTION 0: NOTES AND PRELIMINARIES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
got_labels=0;
Lchol=[];
nochol= 0;
if nargin==6 && ~isempty(labels) 
    got_labels=1;
end

%PARFOR
%RUN FUNCTION THAT TESTS COLLINEARITY?
%PRINT NAIVE SEs for benchmark?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECTION 1: CONSTRUCT DESIGN MATRIX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,time_dummies]                 = unique(year,'rows','stable');
time_dummies 	                   = sparse((1:size(y,1))',time_dummies',1,size(y,1),max(time_dummies));

[~,~,cluster_dummies]              = unique(clusterID,'rows','stable');
cluster_dummies 	               = sparse((1:size(y,1))',cluster_dummies',1,size(y,1),max(clusterID));
X                                  = [D time_dummies cluster_dummies controls];
K                                  = size(X,2);
N_D                                = size(D,2);
N_cluster                          = max(clusterID);
Nobs                               = size(y,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECTION 2: COMPUTE REGRESSION COEFFICIENTS USING FWL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIterations                      = 10000; %iteration for the pcg solver
tol                                = 1e-10; %tol for pcg
xx                                 = X'*X;
xy                                 = X'*y;
try
Lchol 				               = ichol(xx,struct('type','ict','droptol',1e-10,'diagcomp',0.1));
[coeff, flag]		               = pcg(xx,xy,tol,numIterations,Lchol,Lchol');
catch
nochol                             = 1;
[coeff, flag]		               = pcg(xx,xy,tol,numIterations);
end
adj_REGHDFE                        = (Nobs-K-1)/(Nobs-N_cluster-(K-1)); %(N_cluster/(N_cluster-1))*(Nobs-1)/(Nobs-K) %
res                                = sqrt(adj_REGHDFE)*(y-X*coeff); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECTION 3: LEAVE OUT RESIDUALS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meat                               = zeros(K,K);
meat_naive                         = zeros(K,K);
if N_D>100 || Nobs > 10000 %large data, run with parfor
parfor c=1:N_cluster
    
    %% Part 1: Estimate Regression Coefficient After Leaving Cluster c Out
    sel                            = clusterID ~= c;
    y_use                          = y(sel);  
    X_use                          = X(sel,:);
    xx                             = X_use'*X_use;
    xy                             = X_use'*y_use;
    if nochol == 0
    [coeff_leave_out, flag]		   = pcg(xx,xy,tol,numIterations,Lchol,Lchol'); %coefficient vector after leaving a cluster out.
    end
    if nochol == 1
    [coeff_leave_out, flag]		   = pcg(xx,xy,tol,numIterations); %coefficient vector after leaving a cluster out.
    end
    %coeff_leave_out(1:N_D)
    
    %% Part 2: Form KSS sigma_i, compute meat of the sandwich
    sel                            = clusterID==c;
    y_cluster                      = y(sel);
    X_cluster                      = X(sel,:);
    eta_h                          = y_cluster-X_cluster*coeff_leave_out; %leave cluster c out residuals.
    score_left                     = X_cluster'*(y_cluster-mean(y));
    score_right                    = X_cluster'*eta_h;
    meat                           = meat+score_left*score_right'; %X'*sigma_KSS*X for a given cluster
    
    res_cluster                    = res(sel);
    score_left                     = X_cluster'*res_cluster;
    score_right                    = X_cluster'*res_cluster;
    meat_naive                     = meat_naive+score_left*score_right'; %X'*sigma_naive*X for a given cluster
end
else
for c=1:N_cluster
    
    %% Part 1: Estimate Regression Coefficient After Leaving Cluster c Out
    sel                            = clusterID ~= c;
    y_use                          = y(sel);  
    X_use                          = X(sel,:);
    xx                             = X_use'*X_use;
    xy                             = X_use'*y_use;
    if nochol == 0
    [coeff_leave_out, flag]		   = pcg(xx,xy,tol,numIterations,Lchol,Lchol'); %coefficient vector after leaving a cluster out.
    end
    if nochol == 1
    [coeff_leave_out, flag]		   = pcg(xx,xy,tol,numIterations); %coefficient vector after leaving a cluster out.
    end
    %coeff_leave_out(1:N_D)
    
    %% Part 2: Form KSS sigma_i, compute meat of the sandwich
    sel                            = clusterID==c;
    y_cluster                      = y(sel);
    X_cluster                      = X(sel,:);
    eta_h                          = y_cluster-X_cluster*coeff_leave_out; %leave cluster c out residuals.
    score_left                     = X_cluster'*(y_cluster-mean(y));
    score_right                    = X_cluster'*eta_h;
    meat                           = meat+score_left*score_right'; %X'*sigma_KSS*X for a given cluster
    
    res_cluster                    = res(sel);
    score_left                     = X_cluster'*res_cluster;
    score_right                    = X_cluster'*res_cluster;
    meat_naive                     = meat_naive+score_left*score_right'; %X'*sigma_naive*X for a given cluster
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECTION 4: PRINT STANDARD ERRORS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE                                 = zeros(N_D,1);
SE_naive                           = zeros(N_D,1);
xx                                 = X'*X;
for q=1:N_D
    v                              = sparse(q,1,1,K,1);
    if nochol == 0
    [right flag]                   = pcg(xx,v,tol,numIterations,Lchol,Lchol');
    end
    if nochol == 1
    [right flag]                   = pcg(xx,v,tol,numIterations);
    end
    left                           = right';    
    V                              = left*meat*right;
    SE(q)                          = sqrt(V);
    V                              = left*meat_naive*right;
    SE_naive(q)                    = sqrt(V);

end
if got_labels == 1 
    for q=1:N_D 
    tell_me=labels{q};    
    s=['Coefficient on ' tell_me ':  ' num2str(coeff(q))];
    disp(s)
    %s=['Robust "White" Standard Error: ' num2str((SE_naive(q)))];
    %disp(s)
    s=['Leave-Out Standard Error: ' num2str((SE(q)))];
    disp(s)
    s=['******************************************'];
    disp(s);
    end
end
coeff = coeff(1:N_D);


