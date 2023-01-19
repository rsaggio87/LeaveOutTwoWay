%% paths
clc
clear
restoredefaultpath
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay'
path(path,'codes'); %this contains the main LeaveOut Routines.

%% parallel envinr (you can uncommment all of this)
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab

%% read data
namesrc='data/unions.csv'; %path to data
data=importdata(namesrc); %import data
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers
year=data(:,3);
y=data(:,5); % outcome variable
union_status = data(:,4); %0-1 variable on whether a given job-year obs is unionized or not
clear data

%% create year effects 
[~,~,controls] = unique(year,'rows','stable');
controls 	   = sparse((1:size(y,1))',controls',1,size(y,1),max(controls));			
controls(:,end-1)=[]; %omit one year
controls       = [controls union_status];

%% Run variance decomposition for unions non-unions jobs
type_of_algorithm='JLA'; %run random projection algorithm;
[sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_KSS(y,id,firmid,controls,[],type_of_algorithm,250,[],[],[],'results_non_unions',union_status,0); %%calling this will print the variance decomposition for non-union jobs.
[sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_KSS(y,id,firmid,controls,[],type_of_algorithm,250,[],[],[],'results_unions',union_status,1); %%calling this will print the variance decomposition for union jobs.
