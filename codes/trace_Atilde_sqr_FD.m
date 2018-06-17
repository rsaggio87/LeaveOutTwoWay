function trace = trace_Atilde_sqr_FD(Fdelta,F,L,pfun_,NSIM)
%Compute the sum of squared eigenvalues associated to Atilde via
%Hutchinson trick to estimate the trace of large matrices.
if nargin<=4
NSIM=100;
end

%Preallocate
trace=zeros(NSIM,1);



%Simulate trace(B^2) via Hutchinson where B=X*S_xx^(-1)*A*S_xx^(-1)*X'
parfor s=1:NSIM
%inversion    
xsimul=rand(size(Fdelta,1),1);
xsimul=1.*(xsimul>0.5)-1.*(xsimul<=0.5);
[aux flag]=pcg(L,Fdelta'*xsimul,1e-10,1000,pfun_);
aux=F*aux;
A_b=F'*(aux-mean(aux));
[aux, flag]=pcg(L,A_b,1e-10,1000,pfun_);
aux=Fdelta*aux; %%this is B*Y
trace(s)=aux'*aux; 
end

%Results
trace=mean(trace);
end

