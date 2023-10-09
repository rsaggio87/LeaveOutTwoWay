function [coeff, res_F]            = areg_KSS(xy,F,cluster_effects) %%%computes OLS coefficient on F after partialling out cluster fixed effects.

%% Settings for pcg
numIterations                      = 10000; %iteration for the pcg solver
tol                                = 1e-10; %tol for pcg

%% Run using FWL
N_cluster                          = size(cluster_effects,2);
K                                  = size(F,2)+ N_cluster;
d 					               = 1./sum(cluster_effects,1)';
coeff_DF                           = sparse(1:size(d),1:size(d),d) * (cluster_effects'*F);
res_F				               = F-cluster_effects*coeff_DF;
xx					               = res_F'*res_F;
Lchol 				               = ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
[coeff, flag]		               = pcg(xx,xy(N_cluster+1:K) - coeff_DF'*xy(1:N_cluster),tol,numIterations,Lchol,Lchol');
end