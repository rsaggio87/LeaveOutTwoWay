%% paths
clc
clear
restoredefaultpath
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay' %% update this to path where you saved the Github package
path(path,'codes'); %this contains the main LeaveOut Routines.

%% parallel envinr (you can uncommment all of this)
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab

%% some parameters for simulated data
wins_value                      = 15; %event studies binned at -wins_value and wins_value
N_clusters                      = 30; %number of clusters.
N_X                             = 10; %number of extra controls.
T                               = 20; %number of time periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
S                               = 1000;
coeff_sim                       = zeros(S,wins_value*2);
V_sim                           = zeros(S,wins_value*2);
V_sim_naive                     = zeros(S,wins_value*2);
parfor s = 1:S
[y, D, clusterID, year, X ]     = create_DGP(wins_value,N_clusters,N_X,T);
[coeff_sim(s,:),SE,SE_naive]    = KSS_SE(y,D,clusterID,year,X);
V_sim(s,:)                      = SE.^2;
V_sim_naive(s,:)                = SE_naive.^2;
end
C                               = cov(coeff_sim);
disp('True Squared Standard Error')
true_V                          = diag(C)
disp('Average of Squared Standard Error using KSS estimator')
mean(V_sim)'
disp('Average of Squared Standard Error using White SE')
mean(V_sim_naive)'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %% RUN ONE DRAW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y, D, clusterID, year, X ]     = create_DGP(wins_value,N_clusters,N_X,T);
labels                          = cell(10,1);
for k = -wins_value:-2
    labels{k + 6}               = ['Event-Study Coefficient at ' num2str(k)]; % Adjust the index
end

for k = 0:wins_value
    labels{k + 5}               = ['Event-Study Coefficient at ' num2str(k)]; % Adjust the index
end
[coeff, SE]                     = KSS_SE(y,D,clusterID,year,X,labels);
out                             = [coeff,SE]; 
s                               = ['results_AFP.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %FINAL GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run do_file "event_study_with_KSS" to print an event-study graph with
% KSS and "standard" standard errors.

