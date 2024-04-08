%% paths
clc
clear
restoredefaultpath
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay' %% update this to path where you saved the Github package
path(path,'codes'); %this contains the main LeaveOut Routines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %Variance Decomposition within each group type 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
namesrc='data/coll_data_testing.csv'; %path to data (collapsed at the match (worker/firm) level, must be sorted by worker identifier). 
data=importdata(namesrc); %import data. 
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers
y=data(:,3); % outcome variable
group_indicator = data(:,4); %indicator for a group type. this grouping variable must be indexed from 1...G. the code will print variance decomposition for each specific value taken by this variable
Tspell=data(:,5); %number of periods for which a given match is observed in the microdata

clear data

%% options for KSS
type_of_algorithm='JLA'; 
leave_out_level='obs';
simulations_JLA=100;
call_results='results__KSS_collapsed_groups';

%% print KSS-corrected variance decomposition for group g after fitting AKM on the **full** sample
for group_index=1:max(group_indicator)
[tabella(:,group_index)]=leave_out_KSS_coll(y,id,firmid,Tspell,leave_out_level,type_of_algorithm,simulations_JLA,call_results,group_indicator,group_index); %%calling this will print the variance decomposition for non-union jobs.
end


%for each subsample (indexed by the column), tabella reports

% # of Workers in a given subsample
% # of Firms in a given subsample
% # of Person-year obs in a given subsample


% mean outcome in a given subsample
% var outcome in a given subsample


% kss-adjusted std of firm effects in a given subsample
% kss-adjusted std of worker effects in a given subsample
% kss-adjusted correlation of worker,firm effects in a given subsample


% variance of firm effects relative to (residualized) variance outcome.
% variance of person effects relative to (residualized) variance outcome.
% 2x covariance of person,firm effects relative to (residualized) variance outcome.

%residualized variance of the outcome means variance of the outcome after
%taking out the variables in controls (in this example year effects)

%% Code to run "standard" KSS analysis on person-year data /// movers only
if 0 == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %Variance Decomposition within each group type 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
namesrc='data/micro_data_testing_no_stayers.csv'; %path to data. 
data=importdata(namesrc); %import data. these data must be already sorted by worker and year. 
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers
year=data(:,3); %time identifier
y=data(:,4); % outcome variable
group_indicator = data(:,5); %indicator for a group type. this grouping variable must be indexed from 1...G. the code will print variance decomposition for each specific value taken by this variable
clear data

%% options for KSS
type_of_algorithm='JLA'; 
leave_out_level='matches';
simulations_JLA=100;
call_results='results__KSS_groups';

%% print KSS-corrected variance decomposition for group g after fitting AKM on the full-sample
for group_index=1:max(group_indicator)
[tabella(:,group_index)]=leave_out_KSS(y,id,firmid,[],leave_out_level,type_of_algorithm,simulations_JLA,call_results,group_indicator,group_index); %%calling this will print the variance decomposition for non-union jobs.
end
end

%for each subsample (indexed by the column), tabella reports

% # of Workers in a given subsample
% # of Firms in a given subsample
% # of Person-year obs in a given subsample


% mean outcome in a given subsample
% var outcome in a given subsample


% kss-adjusted std of firm effects in a given subsample
% kss-adjusted std of worker effects in a given subsample
% kss-adjusted correlation of worker,firm effects in a given subsample


% variance of firm effects relative to (residualized) variance outcome.
% variance of person effects relative to (residualized) variance outcome.
% 2x covariance of person,firm effects relative to (residualized) variance outcome.

%residualized variance of the outcome means variance of the outcome after
%taking out the variables in controls (in this example year effects)



%% Code to run "standard" KSS analysis on person-year data
if 0 == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %Variance Decomposition within each group type 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data
namesrc='data/kss_subsample_testing_data.csv'; %path to data. 
data=importdata(namesrc); %import data. these data must be already sorted by worker and year. 
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers
year=data(:,3); %time identifier
y=data(:,4); % outcome variable
group_indicator = data(:,5); %indicator for a group type. this grouping variable must be indexed from 1...G. the code will print variance decomposition for each specific value taken by this variable
clear data

%% create year fixed effects
[~,~,controls] = unique(year,'rows','stable');
controls 	   = sparse((1:size(y,1))',controls',1,size(y,1),max(controls));			
controls(:,end-1)=[]; %omit one year

%% options for KSS
type_of_algorithm='JLA'; 
leave_out_level='matches';
simulations_JLA=100;
call_results='results__KSS_groups';

%% print KSS-corrected variance decomposition for group g after fitting AKM on the full-sample
for group_index=1:max(group_indicator)
[tabella(:,group_index)]=leave_out_KSS(y,id,firmid,controls,leave_out_level,type_of_algorithm,simulations_JLA,call_results,group_indicator,group_index); %%calling this will print the variance decomposition for non-union jobs.
end
end

%for each subsample (indexed by the column), tabella reports

% # of Workers in a given subsample
% # of Firms in a given subsample
% # of Person-year obs in a given subsample


% mean outcome in a given subsample
% var outcome in a given subsample


% kss-adjusted std of firm effects in a given subsample
% kss-adjusted std of worker effects in a given subsample
% kss-adjusted correlation of worker,firm effects in a given subsample


% variance of firm effects relative to (residualized) variance outcome.
% variance of person effects relative to (residualized) variance outcome.
% 2x covariance of person,firm effects relative to (residualized) variance outcome.

%residualized variance of the outcome means variance of the outcome after
%taking out the variables in controls (in this example year effects)



