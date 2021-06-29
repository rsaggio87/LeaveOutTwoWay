function [theta,theta_KSS] = kss_quadratic_form(sigma_i,A_1,A_2,beta,Bii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%THIS IS OUR ESTIMATOR, SIMPLE AS THAT!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
right                               = A_2*beta;
left                                = A_1*beta;
theta                               = cov(left,right);
theta                               = theta(1,2);
dof                                 = size(left,1)-1
theta_KSS                           = theta-(1/dof)*sum(Bii.*sigma_i);
end

