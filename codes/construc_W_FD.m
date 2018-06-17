function [W_to_use, my_first_part] = construc_W_FD(ydelta,Fdelta,L,pfun_,A_b,Lambda_B,I_Lambda_P,L_P,eta_h)
%% Constructs C*Y=B*Y-0.5(Lambda_B*eta_h+xi_hat) where xi_hat is the residual when regressing X onto (1-Lambda_p)^(-1)*Lambda_B*Y. 
%% Notice that B*Y=X(X'X)^(-1)A*bOLS where Bols is the OLS estimate. A*bOLS is provided as input
%% This is used to compute the SE of the associated quadratic form

%First Step: B*Y=X(X'X)^(-1)A*bOLS
[aux, flag]=pcg(L,A_b,1e-10,1000,pfun_);
my_first_part=Fdelta*aux; %%this is B*Y

%Second Step: 0.5(Lambda_B*eta_h+xi_hat)
xy=Lambda_B*ydelta;
[ydelta_xi, flag]=pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
xy=Fdelta'*ydelta_xi;
[b, flag]=pcg(L,xy,1e-10,1000,pfun_);
xi=ydelta_xi-Fdelta*b;
my_second_part=0.5*(Lambda_B*eta_h+xi); %% this is 0.5(Lambda_B*eta_h+xi_hat)

%Finish
W_to_use=my_first_part-my_second_part; %%recall the minus!

end

