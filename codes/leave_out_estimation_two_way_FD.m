function [theta, V, sigma_predict, COV_R1, gamma_sq,Fstatistic,b_1,theta_1]= leave_out_estimation_two_way_FD(ydelta,Fdelta,F,L,pfun_,bias_part,Lambda_B,Lambda_P,I_Lambda_P,L_P,eta_h,W_to_use,my_first_part,x1bar,lambda_1)
% % Perform Leave One Out Estimation of Variance Components in AKM in first
%differences.

%% Set up dimensions
Ndelta=size(ydelta,1);
NT=sum(diag(F'*F));
NNarg=13;

if nargout>3 && nargin > NNarg
DO_R1=1;  
aux=lambda_1*(x1bar.^2);
aux=spdiags(aux,0,size(eta_h,1),size(eta_h,1));
Lambda_B2=Lambda_B-aux;
end
if nargout<=3 || nargin <= NNarg
x1bar=0;
lambda_1=0;
DO_R1=0;
Lambda_B2=sparse(size(eta_h,1),size(eta_h,1));    
end
   


%% %% STEP 1: Leave out estimate of the quadratic form
theta=bias_part-(1/NT)*ydelta'*Lambda_B*eta_h;

%% Step 2: Construct SE in the high rank case
%Here we focus on estimation of the sampling variance of theta in the high
%rank case, see end of page 25. We first focus on computing the right part
%of the formula. That part corresponds to the variance of a quadratic form
%with a vector of normal random variables each with variance 
%\tilde{\sigma}^{2}_{i}. We therefore approximate this object via simulations 
%We also proceed in a similar way to compute the object needed to compute
%inference in the case in which q=1.
if nargout>=2
Pii=diag(Lambda_P);
Bii=diag(Lambda_B);
sigma_i=ydelta.*eta_h;
movers=(Pii>0);
stayers=Pii==0;
Pii=Pii(movers,:);
Bii=Bii(movers,:);
bigY=sigma_i(movers);
Nuse=size(Pii,1);
[Xuse, ~, index_unique]=unique([Pii Bii],'rows'); %construct unique combinations
sigma_use = accumarray(index_unique,bigY,[],@mean); %average value of sigma_i within [Bii Pii] cell.
weight=accumarray(index_unique,1);
%% Step 2B: bandwidth
hBest=1/Nuse^(1/3); %choice of constant is arbitrary. We tried several alternatives (e.g. cross validation) and our final estimate seems to change little for different choices of h.  
%% Step 2C: Perform non-parametric fit and provide some diagnostics.
f = fit(Xuse,sigma_use,'lowess','Normalize','on','Weights',weight,'Span',hBest);
sigma_np = feval(f,Xuse);
selnan=isnan(sigma_np);
sigma_np(selnan)=sigma_use(selnan); %%just take the average of sigma_i within cell for cases where we are unable to provide fit
%% Step 2D: Assign
sigma_np=sigma_np(index_unique);
sigma_predict=zeros(size(Fdelta,1),1);
sigma_predict(movers)=sigma_np;
sigma_predict(stayers)=sigma_i(stayers);
%% Step 2E: Construct Right End Part of the matrix (via simulations)
NSIM=round(10e6/Ndelta);
aux_SIM=zeros(NSIM,1);
if nargout>2 && nargin > NNarg
aux_SIM2=zeros(NSIM,1);
end
vsim=randn(Ndelta,NSIM).*movers.*sqrt(sigma_predict); %this will only work on newest versions of matlab. 
parfor s=1:NSIM
       v=vsim(:,s);
       aux=Fdelta'*v;
       [coeff, flag]=pcg(L,aux,1e-5,1000,pfun_); 
       aux=v-Fdelta*coeff;
       bias_part=var(F*coeff)*(NT-1); % This is v'Bv, premultiply by NT to avoid double multiplication wrt last line.
       [aux, flag]=pcg(I_Lambda_P,aux,1e-5,1000,L_P,L_P');
       aux_SIM(s)=bias_part-v'*Lambda_B*aux;
       %Extra bit if we also want calculations for the case when q=1
       if DO_R1==1
       aux_SIM2(s) = bias_part - v'*Lambda_B2*aux - lambda_1*(v'*x1bar)*(x1bar'*v);
       end 
end
left_part=var(aux_SIM);

%% Step 2F: Construct Left End Part of the matrix (no clustering yet)
inner=(W_to_use).*(W_to_use).*sigma_predict;
SUM=sum(inner);
%Complete
V=(1/NT^2)*(4*SUM-left_part);

%% Step 3: Elements do conducts inference in q=1 case
if DO_R1==1
%First thing I need is to form C_2*y=B_2*y-0.5(Lambda_B2*eta_h+xi_hat_2)
%where:  
%xi_hat_2 is the residual when regressing X onto
%(1-Lambda_p)^(-1)*Lambda_B2*y

%B_2=B-lambda_1*x1bar*x1bar'

%Lambda_B2=Lambda_B-diag(lambda_1*x1bar*x1bar') --> defined at the
%beginning.



%start with the first piece: 
%B_2*y=(B-lambda_1*x1bar*x1bar')y

%Recall that B*y=my_first_part
b_1=sum(x1bar.*ydelta);
first_piece=my_first_part-(lambda_1*x1bar)*b_1;


%Now the second bit
xy=Lambda_B2*ydelta;
[ydelta_xi, flag]=pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
xy=Fdelta'*ydelta_xi;
[b, flag]=pcg(L,xy,1e-10,1000,pfun_);
%                 if flag > 0
%                     error(['PCG FLAG: ' num2str(flag)])
%                 end
xi=ydelta_xi-Fdelta*b;
second_piece=0.5*(Lambda_B2*eta_h+xi);
W_2=first_piece-second_piece; %this gives me W_2(i)=\sum_{l \neq i}C_{il}y_{l}

%Let's set-up things to calculate the VCM
COV_R1=zeros(2,2);

%Start with COV(1,1)
%COV_R1(1,1)=sum(x1bar.^2.*sigma_i); %do not divide
COV_R1(1,1)=sum(x1bar.^2.*sigma_predict); %do not divide

%Next is COV(1,2)
COV_R1(1,2)=2*sum(x1bar.*sigma_predict.*W_2)/(NT);
COV_R1(2,1)=COV_R1(1,2);

%More complicated step is COV(2,2). Similar structure of calculations as
%performed when q=0 --> variance of a quadratic form with mean zero normal
%vector. See the calculations computed there (inefficient to write another
%parfor).

second_bit=var(aux_SIM2);
COV_R1(2,2)=(4*sum(W_2.*W_2.*sigma_predict)-second_bit)/(NT^2);

%Gamma Squared
gamma_sq=((lambda_1^2/NT^2)*(COV_R1(1,1)^2))/(COV_R1(2,2));

%Fstatistic
Fstatistic=(b_1^2)/COV_R1(1,1);

%For AM Method
theta_1=theta-(lambda_1/(NT))*(b_1^2-COV_R1(1,1));
end



end
end

    


