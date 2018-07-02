function [trace_fe,trace_cov,trace_pe] = trace_Atilde_sqr(X,F,D,xx,Lchol,NSIM)
%Compute the sum of squared eigenvalues associated to Atilde via
%Hutchinson trick to estimate the trace of large matrices.

%Read
do_cov=1;
do_pe=1;

if nargout<=1
    do_cov=0;
    do_pe=0;
end

if nargout<=2
    do_pe=0;
end

if nargin<=5
NSIM=100;
end

N=size(D,2);
J=size(F,2); %F should be F=F*S;
K=size(xx,2)-N-J;

%Preallocate
trace_fe=zeros(NSIM,1);
trace_cov=zeros(NSIM,1);
trace_pe=zeros(NSIM,1);

parfor s=1:NSIM

%inversion    
xsimul=rand(size(X,1),1);
xsimul=1.*(xsimul>0.5)-1.*(xsimul<=0.5);
[coeff flag]=pcg(xx,X'*xsimul,1e-10,1000,Lchol,Lchol');


%Variance of firm effects
eff=F*coeff(N+1:N+J,1);
if K>0
A_b=[zeros(N,1); F'*(eff-mean(eff)); zeros(K,1)];
end
if K==0
A_b=[zeros(N,1); F'*(eff-mean(eff))];
end
[aux, flag]=pcg(xx,A_b,1e-10,1000,Lchol,Lchol');
aux=X*aux;
trace_fe(s)=aux'*aux;

%Variance of person effects
if do_pe==1
eff=D*coeff(1:N,1);
if K>0
A_b=[D'*(eff-mean(eff)); sparse(J,1); zeros(K,1)];
end
if K==0
A_b=[D'*(eff-mean(eff)); sparse(J,1)];
end
[aux, flag]=pcg(xx,A_b,1e-10,1000,Lchol,Lchol');
aux=X*aux;
trace_pe(s)=aux'*aux;
end
%Covariance of person,firm effects
if do_cov==1
pe=D*coeff(1:N,1);
fe=F*coeff(N+1:N+J,1);
if K>0
A_b=[0.5*D'*(fe-mean(fe)); 0.5*F'*(pe-mean(pe)); zeros(K,1)];
end
if K==0
A_b=[0.5*D'*(fe-mean(fe)); 0.5*F'*(pe-mean(pe))];
end
[aux, flag]=pcg(xx,A_b,1e-10,1000,Lchol,Lchol');
aux=X*aux;
trace_cov(s)=aux'*aux;
end
end

%Results
trace_fe=mean(trace_fe);
trace_cov=mean(trace_cov);
trace_pe=mean(trace_pe);
end

