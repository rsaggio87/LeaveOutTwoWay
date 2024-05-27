%% paths
clc
clear
restoredefaultpath
cd '/Users/raffaelesaggio/Documents/GitHub/LeaveOutTwoWay' %% update this to path where you saved the Github package
path(path,'codes'); %this contains the main LeaveOut Routines.

%% parallel envinr (you can uncommment all of this)
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %TABLE 2: Variance Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
%% read data
namesrc='data/unions_SAKM.csv'; %path to data
data=importdata(namesrc); %import data
id=data(:,1); %worker identifiers (unrestricted interaction of worker id and union dummy)
firmid=data(:,2); %firm identifiers 
id_orig=data(:,3); %original worker id identifier
year=data(:,4);
y=data(:,5); % outcome variable
union_status = data(:,6); %0 for non-union job; 1 for union job.
age=data(:,7);
clear data

%% age polynomials
age2 = ((age-40)/40).^2;
age3 = ((age-40)/40).^3;

%% create year effects 
[~,~,controls] = unique(year,'rows','stable');
controls 	   = sparse((1:size(y,1))',controls',1,size(y,1),max(controls));			
controls(:,end-1)=[]; %omit one year
controls       = [controls age2 age3];

%% options for KSS
type_of_algorithm='JLA'; 
leave_out_level='matches';
simulations_JLA=250;
call_results='SAKM_RESULTS_LEAVE_OUT_CS';

%% compute variance decomposition 
[cov_pe var_pe_union var_pe_non_union]= leave_out_KSS_SAKM(y,id,firmid,id_orig,controls,leave_out_level,type_of_algorithm,simulations_JLA,call_results,union_status,year); %%calling this will print the variance decomposition for non-union jobs.
end


