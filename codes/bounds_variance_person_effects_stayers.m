function [sigma2_alpha_LB, sigma2_alpha_UB] = bounds_variance_person_effects_stayers(id,X_pe,Bii_pe,b,y,sigma_i);
T                    = accumarray(id,1);
stayers              = T==1;
stayers              = stayers(id);
[~,  sigma2_alpha_UB]= kss_quadratic_form(sigma_i,X_pe,X_pe,b,Bii_pe);
MSE                  = (y-mean(y)).^2; 
sigma_i(stayers)     = MSE(stayers);
[~, sigma2_alpha_LB] = kss_quadratic_form(sigma_i,X_pe,X_pe,b,Bii_pe);
end

