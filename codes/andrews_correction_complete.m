function [var_corrected_fe, var_corrected_pe,var_corrected_cov] = andrews_correction_complete(y,F,D,controls,NSIM)
%% Author: Raffaele Saggio.
%% email:  raffaele.saggio@berkeley.edu

%% Description
%This code performs Andrews et al. (2008) homoskedastic correction in the 
%context of a two way fixed effects model a la AKM.
%The input data should correspond to a connected set to insure proper
%identification.
%To compute the trace needed to perform the correction we follow 
%Gaure (2014) and in particular implement the Hutchinson trick. 

%% Description of Inputs
%y: (NTx1) vector of outcomes, where NT=total number of person year obs.
%D: (NTxN) matrix of workers' assignments.
%F: (NTxJ) matrix of firms' assignements.
%Controls: (NTxK) matrix of controls.
%NSIM: number of simulations to perform hutchinson trick. Default is 100.

%% Description of Output
%var_corrected_fe: corrected variance of firm effects (person year
%weighted)
%var_corrected_pe: corrected variance of worker effects (person year
%weighted)
%corrected covariance of worker, firm effects (person year
%weighted)

%% Additional notes:
%This function is intended for general use on two-way fixed effects model.
%If the user already calculated (Bii) for all i for a given variance
%decomposition parameter (e.g. variance of firm effects) then one can use
%knowledge of this term to compute the Andrews correction without the need
%to run the simulation to compuate an estimate of the trace. 



no_controls=0;
if nargin==3
controls=[]; 
no_controls=1;
NSIM=100;    
end    

if nargin==4  
NSIM=100;    
end 

%% Step 0: Set Up dimensions
NT=size(y,1);
J=size(F,2);
N=size(D,2);
if no_controls==0
K=size(controls,2);
end
if no_controls==1
K=0;
end
S=speye(J-1);
S=[S;sparse(-zeros(1,J-1))];  %NT+JxNT+J-1 restriction matrix 
%% Step 1: Perform AKM
if no_controls==0
X=[D,F*S,controls];
end
if no_controls==1
X=[D,F*S];
end
xx=X'*X;
xy=X'*y;
L=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
b=pcg(xx,xy,1e-10,1000,L,L');
ahat=b(1:N);
ghat=b(N+1:N+J-1);
xb=X*b;
r=y-xb;
dof=NT-size(xx,1)-1;
sigma_mse=(sum(r.^2)/dof);
%% Step 2: Correct
Vsel_pe=zeros(NSIM,1);
Vsel_fe=zeros(NSIM,1);
Vsel_cov=zeros(NSIM,1);
parfor s=1:NSIM
xsimul=rand(size(xx,1),1);
x=1.*(xsimul>0.5)-1.*(xsimul<=0.5); %Rademacher weights
[first flag]=pcg(xx,x,1e-8,1000,L,L');
%% now keep only the relevant solution from pcg and from drawn x's: Var(FE)
first_fe=first(N+1:N+J-1,:);
vecx_fe=x(N+1:N+J-1);
second_fe=F*S*vecx_fe-mean(F*S*vecx_fe);
COV=cov(second_fe,F*S*first_fe);
Vsel_fe(s)=COV(1,2);
%% now keep only the relevant solution from pcg and from drawn x's: Var(PE)
first_pe=first(1:N,:);
vecx_pe=x(1:N);
second_pe=D*vecx_pe-mean(D*vecx_pe);
COV=cov(second_pe,D*first_pe);
Vsel_pe(s)=COV(1,2);
%% now keep only the relevant solution from pcg and from drawn x's: COV(PE,FE)
first_cov=first(N+1:N+J-1,:);
vecx_cov=x(1:N);
second_cov=D*vecx_cov-mean(D*vecx_cov);
COV=cov(second_cov,F*S*first_cov);
Vsel_cov(s)=COV(1,2);
end
%% Don't foreget to multiply by MSE
correction_fe=sigma_mse*(mean(Vsel_fe));
correction_pe=sigma_mse*(mean(Vsel_pe));
correction_cov=sigma_mse*(mean(Vsel_cov));
%% Get uncorrected variance or covariance.
var_uncorrected_fe=var(F*S*ghat);
var_uncorrected_pe=var(D*ahat);
var_uncorrected_cov=cov(D*ahat,F*S*ghat);
var_uncorrected_cov=var_uncorrected_cov(1,2);
%% calculate adjusted variance.
var_corrected_fe=var_uncorrected_fe-correction_fe;
var_corrected_pe=var_uncorrected_pe-correction_pe;
var_corrected_cov=var_uncorrected_cov-correction_cov;
end