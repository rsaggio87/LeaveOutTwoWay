function [UB,LB,C] = AM_CI(NT,lambda_1,gamma_sq,COV_R1,bar_Beta_1,theta_2)
%Redefine as in pdf
lambda_1=lambda_1/NT;
gamma=sqrt(gamma_sq);
sigma_sq=COV_R1(1,1);
tau_sq=COV_R1(2,2);
rho=COV_R1(1,2)/(sqrt(sigma_sq)*(sqrt(tau_sq)));
C=2*abs(gamma)/(sqrt(1-rho^2));


%Step 1: search for corresponding quantile
load('tabulation_10K')
Cgrid=tabula(:,1);
quant_simul=tabula(:,2);
dist=abs(C-Cgrid);
sel=(dist==min(dist));
z_crit=quant_simul(sel);

%Step 2: find roots of crazy polynomial.
p1=(4*gamma_sq)/(sigma_sq);
p2=((4*gamma*rho)/(sqrt(sigma_sq)))-(8*bar_Beta_1*gamma_sq/(sigma_sq));
p3=(1 -4*gamma_sq*z_crit^2 +(4*gamma_sq*bar_Beta_1^2)/(sigma_sq) - (8*bar_Beta_1*rho*gamma)/(sqrt(sigma_sq)));
p4=(4*bar_Beta_1^2*gamma*rho)/(sqrt(sigma_sq))-2*bar_Beta_1-4*gamma*z_crit^2*rho*sqrt(sigma_sq);
p5=bar_Beta_1^2-rho^2*sigma_sq*z_crit^2;
p=[p1 p2 p3 p4 p5];
b=roots(p);
b=b(imag(b)==0);

%Step 3: conclude
den=(1-rho^2+((2*gamma*b)/(sqrt(sigma_sq))+rho).^2).^(1/2);
LB_s=lambda_1*b.^2+theta_2-rho*sqrt(tau_sq/sigma_sq)*(bar_Beta_1-b)-(z_crit*sqrt(tau_sq)*(1-rho^2))./den;
UB_s=lambda_1*b.^2+theta_2-rho*sqrt(tau_sq/sigma_sq)*(bar_Beta_1-b)+(z_crit*sqrt(tau_sq)*(1-rho^2))./den;
LB=min(LB_s);
UB=max(UB_s);

end

