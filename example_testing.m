%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
%This m-file provides five examples that show how to conduct inference on linear
%combination of linear regression parameters. Such inference is robust to
%heteroskedasticity, many regressors asymptotics and serial correlation in the 
%error term - see Proposition 1 and Remark 9 of Kline, Saggio, Solvsten (KSS henceforth).

%All examples below refer to a case where the user in a prior step has 
%performed KSS leave out estimation of a two-way fixed effects model. This 
%is because "lincom_KSS" is essentially a post-estimation command, optimized 
%to run following the command "leave_out_COMPLETE". 
%Nevertheless, "lincom_KSS" can be applied to ANY linear regression model.

%The user should have installed all the files contained in the GitHub
%repository at https://github.com/rsaggio87/LeaveOutTwoWay to properly run
%the code below.

%See the documentation inside "lincom_KSS" for further details. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Testing whether Firm Fixed Effects are different across two regions.

%In this example we will start by loading a person-year file that contains
%wage information from two time periods on two Italian regions. 
%The object is to test whether the firm effects from these two regions are 
%statistically different.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%Paths
path(path,'~/matlab_bgl/'); %note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
path(path,'codes'); %this contains the main LeaveOut Routines.

%Parpooling
try
%pool=parpool(32,'IdleTimeout', Inf);
[a b]=system('sinfo -p low -o "%C"');
cores=str2num([b(19) b(20)]);
cores=min(cores,64);
pool = parpool('dcs', cores);
end

%Log File
placelog='logs/';
namelog='test_example1';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%How to name the results
nameFile='test_example1';
placeFile='results/';
filename=[placeFile nameFile];

%Import data 
namesrc='src/file_example1.csv'; %where original data is
data=importdata(namesrc);
id=data(:,1); 
firmid=data(:,2);
y=data(:,5);
region=data(:,4); %Region indicator. Value -1 for region 1, Value 1 for region 2.

%Run Leave Out Estimation
leave_out_level='obs';
resid_controls=2; 
andrews_estimates=0;
eigen_diagno=0;
subsample_llr_fit=0;
restrict_movers=0;
do_SE=0;
type_of_algorithm='exact';
sigma2_psi = leave_out_COMPLETE(y,id,firmid,leave_out_level,region,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,[],filename);

%Load data associated with Leave One Out Connected Set. 
results=importdata([filename '.csv']);
y=results(:,1);
firmid=results(:,2);
id=results(:,3);
region=results(:,7:end);
clear results

%Load the matrix of statistical leverages saved by leave_out_COMPLETE
Lambda_P=load([filename '_Lambda_P']);
Lambda_P=Lambda_P.Lambda_P;

%Now create the X associated with the two-way model (dropping as usual last firm).
NT=size(y,1);
D=sparse(1:NT,id',1); 
F=sparse(1:NT,firmid',1);
N=size(D,2);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];
F=F*S;
X=[D,F]; %Design Matrix.

%Test linear combination v'beta where beta=(X'*X)X'*y.
Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
region_dummy=region;
region_dummy(region_dummy==-1)=0;
Z=region_dummy; 
labels={'Region 2 Dummy'};

%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
%function "lincom_KSS" always adds a constant to the user specified matrix
%"Z".

%Therefore in this case the second element of v'*beta returns the regression
%coefficient that captures the difference in firm effects between region 2
%and region 1. 

%We want to conduct robust inference on this particular
%regression coefficient. The function lincom_KSS outputs the t-test
%of this particular linear combination.

%Run
t_stat_example_1=lincom_KSS(y,X,Z,Transform,[],Lambda_P,labels); %Note: the [] implies that we will provide heteroskedatic robust inference (i.e. no clustering) which is consistent with the prior step where we conducted leave out estimation by leaving a single observation out each time (leave_out_level='obs';)

%Close
diary off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Example 2: Projecting firm effects on observables

% In this example, we show how to conduct inference on the
% projection of firm fixed effects onto the variables: worker gender, 
% worker age, firm size, firm age.
% 

% The data come from an unbalanced panel with at most 6 time periods per worker.
% As before, we will first start by conducting leave out estimation of the model
% to obtain: y,X, and the statistical leverages belonging to the leave
% out connected set (where we know the leave-out model is identified). 

% Then, we'll apply the function "lincom_KSS".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Log File
placelog='logs/';
namelog='test_example2';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%Name the file to be saved
nameFile='test_example2';
placeFile='results/';
filename=[placeFile nameFile];

%Import data 
namesrc='src/file_example2.csv'; %where original data is
data=importdata(namesrc);
id=data(:,1);
firmid=data(:,2);
y=data(:,4);
extra_vars=data(:,5:end); %Column 1: gender (1 for female). %Column 2: Workers' Age % Column 3: Log Firm Size %Column 4: Firm Age.

%Run Leave Out Estimation 
leave_out_level='obs';
restrict_movers=1; %this is turn on to make the sample from Example 2 equivalent to the one used in Example 3.
[sigma2_psi] = leave_out_COMPLETE(y,id,firmid,leave_out_level,extra_vars,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE,type_of_algorithm,[],filename);

%Load the .csv file created by the function  leave_out_COMPLETE
results=importdata([filename '.csv']);
y=results(:,1);
firmid=results(:,2);
id=results(:,3);
extra_vars=results(:,7:end);
clear results

%Load the matrix Lambda_P saved by leave_out_COMPLETE
Lambda_P=load([filename '_Lambda_P']);
Lambda_P=Lambda_P.Lambda_P;

%Now create the X associated with the two-way model (dropping as usual last firm).
NT=size(y,1);
D=sparse(1:NT,id',1); 
F=sparse(1:NT,firmid',1);
N=size(D,2);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];
F=F*S;
X=[D,F]; %Design Matrix.

Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
Z=extra_vars;
labels={'Gender';'Worker Age'; 'Log Firm Size';'Firm Age'};

%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
%function "lincom_KSS" always adds a constant to the user specified matrix
%"Z". Also recall that beta=(X'*X)^(-1)X'y;

%Therefore in this case v'*beta in Remark 9 of KSS gives back the 4 (gender,
%worker age, log firm size, firm age) t-statistics obtained when projecting
%the firm effects onto a constant, gender dummy, Age, Log Firm Size and Firm Age.

t_stat_example_2=lincom_KSS(y,X,Z,Transform,[],Lambda_P,labels);

%Close
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 3: Projecting firm effects on observables (clustered by "match")

% This example mimics example 2 but allows arbitrary serial correlation of errors
% within a worker-firm match.

% To accommodate clustering, we re-run the command "leave_out_COMPLETE"
% this time leaving out the entire match to conduct estimation.

% Next, we apply "lincom_KSS" making sure to turn on the cluster option.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Log File
placelog='logs/';
namelog='test_example3';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%Name the file to be saved
nameFile='test_example3';
placeFile='results/';
filename=[placeFile nameFile];

%Import data 
namesrc='src/file_example2.csv'; %where original data is
data=importdata(namesrc);
id=data(:,1);
firmid=data(:,2);
y=data(:,4);
extra_vars=data(:,5:end); %Column 1: gender (1 for female). %Column 2: Workers' Age % Column 3: Log Firm Size %Column 4: Firm Age.

%Run Leave Out Estimation 
leave_out_level='matches';
restrict_movers=1 %this is especially important here because for stayers AKM model is not identified when leaving an entire match out!
[sigma2_psi] = leave_out_COMPLETE(y,id,firmid,leave_out_level,extra_vars,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE,type_of_algorithm,[],filename);

%Load the .csv file created by the function  leave_out_COMPLETE
results=importdata([filename '.csv']);
y=results(:,1);
firmid=results(:,2);
id=results(:,3);
extra_vars=results(:,7:end);
[~,~,match_id]=unique([id firmid],'rows');
clear results

%Load the matrix Lambda_P saved by leave_out_COMPLETE (Note: now this matrix is going to be block-diagonal instead of diagonal as in Example 2) 
Lambda_P=load([filename '_Lambda_P']);
Lambda_P=Lambda_P.Lambda_P;

%Now create the X associated with the two-way model (dropping as usual last firm).
NT=size(y,1);
D=sparse(1:NT,id',1); 
F=sparse(1:NT,firmid',1);
N=size(D,2);
J=size(F,2);
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];
F=F*S;
X=[D,F]; %Design Matrix.

Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
Z=extra_vars;
labels={'Gender';'Worker Age'; 'Log Firm Size';'Firm Age'};

%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
%function "lincom_KSS" always adds a constant to the user specified matrix
%"Z". Also recall that beta=(X'*X)^(-1)X'y;

%Therefore in this case v'*beta in Remark 9 of KSS gives back the 4 (gender,
%worker age, log firm size, firm age) t-statistics obtained when projecting
%the firm effects onto a constant, gender dummy, Age, Log Firm Size and Firm Age.

t_stat_example_3=lincom_KSS(y,X,Z,Transform,match_id,Lambda_P,labels); %cluster the standard errors by match identifiers now.

%Close
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 4: Joint Test (cluster, heteroskedatic robust)

% Here we show how to run a joint test in the setup highlighted in example
% 3. In particular, we will report the statistic and associated p-value
% associated with the hypothesis that the variables contained in Z 
% (gender, age, firm size, firm age)jointly predict the firm effects 
% from the AKM model.

% All the relevant matrices are stored in memory, the user needs to run
% this function after running Example 3. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Log File
placelog='logs/';
namelog='test_example4';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%Run
[~,~,~,~,~,stat_4,pvalue_4]=lincom_KSS(y,X,Z,Transform,match_id,Lambda_P,labels);

%Close
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 5: A Different Joint-Test (cluster, heteroskedatic robust)

% Here we provide a slightly different joint test. Instead of testing
% whether all the variables in Z jointly predict the firm effects, which is
% what the function lincom_KSS outputs by default, we test
% whether firm size and firm age are jointly significant in a regression 
% where the firm effects are regressed onto a constant, gender, age, firm
% size and firm age.

% This example shows how, by appropriately specifying the matrix "restrict"
% one can conduct such test.

% All the relevant matrices are stored in memory, the user needs to run
% this function after having completed Example 3. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Log File
placelog='logs/';
namelog='test_example5';
logname=[placelog namelog '.log'];
system(['rm ' logname])
diary(logname)

%Build the sub-restriction
%Remember in the function Z will have the following structure 
%
%            Z = [constant gender age firm_size firm_age]
%
%We want to test whether  worker age and firm_age have jointly predict the 
%firm effects. 
%To do so, we create the appropriate matrix "restrict" so that 
%
%           v= restrict*(Z'Z)^(-1)Z'Transform
%


restrict=[0  0  0  1  0; 0  0  0  0  1];
[~,~,~,~,~,stat_5,pvalue_5]=lincom_KSS(y,X,Z,Transform,match_id,Lambda_P,labels,restrict);


%Close
diary off
delete(pool)