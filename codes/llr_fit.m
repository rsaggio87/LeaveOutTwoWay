function sigma_predict= llr_fit(Lambda_P,Lambda_B,y,eta_h,subsample_llr_fit,sizeRandomSample)
if nargin==5
sizeRandomSample=0.20;    
end
%% Set up dimensions
NT=size(eta_h,1);
%% Run LLR
tic
sigma_i=y.*eta_h;
Pii=diag(Lambda_P);
Bii=diag(Lambda_B);
hBest=1/NT^(1/3);
if subsample_llr_fit == 0 
f = fit([Pii Bii],sigma_i,'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,[Pii, Bii]);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);
end
if subsample_llr_fit == 1
sel=(rand(NT,1)<=sizeRandomSample);    
f = fit([Pii(sel) Bii(sel)],sigma_i(sel),'lowess','Normalize','on','Span',hBest);
sigma_predict = feval(f,[Pii, Bii]);
selnan=isnan(sigma_predict);
sigma_predict(selnan)=sigma_i(selnan);
end
disp('Time to Perform non-parametric fit')
toc
disp('Fraction of negative predicted values for tilde sigma_i^2')
mean(sigma_predict<0)
end    


