%% this is a m-file created to show how one can obtain unweighted KSS
%% adjusted variance of firm effects.


%% setting up
clc
clear
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay'
path(path,'codes'); %this contains the main LeaveOut Routines. path(path,'CMG'); % CMG package http://www.cs.cmu.edu/~jkoutis/cmg.html 
%[result,output] = evalc('installCMG(1)'); %installs CMG routine (silently) 
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores 
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab

%% Import data.
namesrc='data/black_soc3_firm.csv'; %path to original testing data 
data=importdata(namesrc); %import data
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers
y=data(:,3); % outcome variable

clear data
%% Run. Works only for variance of firm effects.
[sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_KSS(y,id,firmid,[],'obs','exact');