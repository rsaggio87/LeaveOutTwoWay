function [Lambda_P, Lambda_B_fe, Lambda_B_cov, Lambda_B_pe] = eff_res(X,xx,Lchol,N,J,elist,clustering_level)
%This function calculates, using parallel coding, (Pii,Bii) for a general 
%two way fixed effects model. The code is likely to be slow on large
%datasets.


%Dimensions
NT=size(X,1);
M=size(elist,1);

%PreCreate
Pii=zeros(M,1);
Bii_fe=zeros(M,1);
Bii_cov=zeros(M,1);
Bii_pe=zeros(M,1);

%Options for solver
numIterations = 300; %iteration for the pcg solver
tol=1e-6; %tol for pcg

%slice
Xright=X(elist(:,2),:);
Xleft=X(elist(:,1),:);

%Loop
if ~strcmp(clustering_level,'obs')
parfor i=1:M     
[xtilde_right, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol,Lchol');
[xtilde_left, flag]= pcg(xx,Xleft(i,:)',tol,numIterations,Lchol,Lchol');

%Statistical Leverage
Pii(i)=Xleft(i,:)*xtilde_right;

%Bii for Variance of Firm Effects
aux_right=xtilde_right(N+1:N+J-1,:);
aux_left=xtilde_left(N+1:N+J-1,:);
COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
Bii_fe(i)=COV(1,2)*(NT-1);

%Bii for Variance of Person Effects
aux_right=xtilde_right(1:N);
aux_left=xtilde_left(1:N);
COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
Bii_pe(i)=COV(1,2)*(NT-1);

%Bii for Covariance of Person, Firm Effects
aux_right=xtilde_right(N+1:N+J-1);
aux_left=xtilde_left(1:N);
COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
Bii_cov(i)=COV(1,2)*(NT-1);
end
end

if strcmp(clustering_level,'obs')
parfor i=1:M     
[xtilde, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol,Lchol');

%Statistical Leverage
Pii(i)=Xright(i,:)*xtilde;

%Bii for Variance of Firm Effects
aux_right=xtilde(N+1:N+J-1,:);
aux_left=xtilde(N+1:N+J-1,:);
COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
Bii_fe(i)=COV(1,2)*(NT-1);

%Bii for Variance of Person Effects
aux_right=xtilde(1:N);
aux_left=xtilde(1:N);
COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
Bii_pe(i)=COV(1,2)*(NT-1);

%Bii for Covariance of Person, Firm Effects
aux_right=xtilde(N+1:N+J-1);
aux_left=xtilde(1:N);
COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
Bii_cov(i)=COV(1,2)*(NT-1);
end
end


%Create the matrices.
rows=elist(:,1);
column=elist(:,2);
%Lambda P
Lambda_P=sparse(rows,column,Pii,NT,NT);
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.
%Lambda B var(fe)
Lambda_B_fe=sparse(rows,column,Bii_fe,NT,NT);
Lambda_B_fe=Lambda_B_fe+triu(Lambda_B_fe,1)'; %make it symmetric.
%Lambda B cov(fe,pe)
Lambda_B_cov=sparse(rows,column,Bii_cov,NT,NT);
Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)';
%Lambda B, var(pe)
Lambda_B_pe=sparse(rows,column,Bii_pe,NT,NT);
Lambda_B_pe=Lambda_B_pe+triu(Lambda_B_pe,1)'; %make it symmetric.
end

