%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This m-file is used to compute leave out estimates on a two-way 
%fixed effects model using a test sample.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %External Packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'~/matlab_bgl/'); 
%note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
%note: MakeCMG called only if running leave_out_FD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %PARALLEL ENVIRONMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%Note: You should decide which set-up is most suitable given
%      your own setup. I made use of the dcs configuration
%      available for the Berkeley hpc. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pool = parpool('dcs', 64);
%pool=parpool(12);
%pool=parpool(str2num(getenv('SLURM_CPUS_PER_TASK')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %CALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%Note: This is where you should specify where is the src .csv file
%      and how you would like to call the log-file created by the
%      leave out functions

%Make sure that the src .csv file is sorted by worker id and year
%(xtset id year).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namesrc='../src/test.csv';
namelog='test';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=importdata(namesrc);
id=data(:,1);
firmid=data(:,2);
year=data(:,3);
y=data(:,4);
controls=data(:,5:end); %in the test data this is a matrix containing year dummies, omitting one year (1995 in this case). Check description of leave_out_FD to understand how we handle these controls.
clear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %RUN LEAVE-OUT (GENERAL BUT SLOW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 == 1
%Options (see description inside leave_out_COMPLETE)
leave_out_level='workers';
andrews_estimates=0;
eigen_diagno=0;
subsample_llr_fit=0;
restrict_movers=0;
resid_controls=1;
Ndiagno=1;
%Log File
logname=['../logs/leave_one_out_complete' namelog];
system(['rm ' logname])
diary(logname)
[sigma2_psi, V] = leave_out_COMPLETE(y,id,firmid,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,Ndiagno,namelog);    
diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %RUN LEAVE-OUT (JUST VARIANCE OF FIRMS EFFECTS BUT FAST)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
%Options (see description inside leave_out_FD)
leave_out_level='workers';
type_algorithm='JL';
eigen_diagno=1;
eigen_fast=0; %This checks the conditions of theorem 1 of KSS. It is computationally intensive so I turned this off but I suggest turning this option to 1 at one point.  
do_montecarlo=0;
%Log File
logname=['../logs/leave_one_out_firm_effects_only' namelog];
system(['rm ' logname])
diary(logname)
[sigma2_psi, V]= leave_out_FD(y,id,firmid,leave_out_level,controls,type_algorithm,eigen_diagno,eigen_fast,do_montecarlo,namelog);    
diary off
end
%close
%delete(pool)
