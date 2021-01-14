function [sigma2_alpha_LB, sigma2_alpha_UB] = bounds_variance_person_effects_stayers(id,peso,Mii,eta,X_pe,Bii_pe,b,outcome);
mean(y)
save('cacca2')
T                   = accumarray(id,1);
stayers             = T==1;
stayers             = stayers(id);
Mii(stayers)        = 1./peso(stayers); %estimate on the variance of person effects will only be heteroskedastic robust
eta_h				= eta./Mii; %Leave one out residual
sigma_i				= y.*eta_h; %KSS estimate of individual variance.
[~,  sigma2_alpha_UB]= kss_quadratic_form(sigma_i,X_pe,X_pe,b,Bii_pe);
[~, sigma2_alpha_LB]= kss_quadratic_form(y.^2,X_pe,X_pe,b,Bii_pe);
end

