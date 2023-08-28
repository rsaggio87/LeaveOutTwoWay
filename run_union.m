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
firmid=data(:,2); %firm identifiers (unrestricted interaction of firmid and union dummy)
firmid_orig=data(:,3); %original firm identifiers
year=data(:,4);
y=data(:,5); % outcome variable
union_status = data(:,6); %0 for non-union job at firm where union share is 0%; 1 for union job; 2 for non-union job at firm with non-zero share of union jobs
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
call_results='results__KSS_unions';

%% run
for type_decomposition=0:2 
tabella(:,type_decomposition+1) = leave_out_KSS(y,id,firmid,firmid_orig,controls,leave_out_level,type_of_algorithm,simulations_JLA,call_results,union_status,type_decomposition); %%calling this will print the variance decomposition for non-union jobs.
end

%% print table into a .csv
tabella=[tabella(:,2) tabella(:,3) tabella(:,1)];
tabella=full(tabella);
s=['tables/VARIANCE_DECOMPOSITION.csv'];
dlmwrite(s, tabella, 'delimiter', '\t', 'precision', 16);