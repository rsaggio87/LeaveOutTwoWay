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
%Note: The user should decide which set-up is most suitable given
%      her own configuration. Make sure to delete the pool once
%      estimation has been carried out. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
pool = parpool('dcs', 64);
%delete(gcp('nocreate'))
%pool=parpool(32);
%pool=parpool(str2num(getenv('SLURM_CPUS_PER_TASK')));
end
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
controls=data(:,5:5); %in the test data this is a matrix containing year dummies, omitting one year (2001 in this case) in order to avoid collinearity.
clear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %RUN LEAVE-OUT (GENERAL BUT SLOW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
leave_out_level='obs';
andrews_estimates=0;
eigen_diagno=0;
restrict_movers=0;
resid_controls=0;
controls=[]; %equivalent to no controls
do_SE=0;
subsample_llr_fit=0;
logname=['../logs/leave_one_out_' leave_out_level '_' namelog '.log'];
system(['rm ' logname])
diary(logname)
sigma2_psi = leave_out_COMPLETE(y,id,firmid,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE,namelog);    
diary off
end
%close
delete(pool)
