function [v,lambda_eig] = eigAux(type_quadratic_form,xx,Lchol,F,D,K)
%This finds the eigenvalues/vector of the matrix Atilde without the need to
%store the matrix Atilde but by looking at the eigenvalues of two
%particular matrices, each of them easy to store and create. 

N=size(D,2);
J=size(F,2); %F should F*S if working with grounded Laplacian
NT=size(D,1);
Sxxdot=[xx sparse(N+J+K,1); sparse(1,N+J+K) 1];


%Variance of Firm Effects
if strcmp(type_quadratic_form,'fe')
    a=[sparse(NT,N) F sparse(NT,K)];
end

%Variance of Person Effects
if strcmp(type_quadratic_form,'pe')
    a=[D sparse(NT,J) sparse(NT,K)];
end

%Covariance of Person,Firm Effects
if strcmp(type_quadratic_form,'cov')
    a=[D sparse(NT,J) sparse(NT,K)];
    b=[sparse(NT,N) F sparse(NT,K)];
end


if ~strcmp(type_quadratic_form,'cov')
abar=sum(a,1)'/(sqrt(NT));
S_aa=a'*a;
ainv=pcg(xx,abar,1e-5,1000,Lchol,Lchol');
Adot=[S_aa -S_aa*ainv;abar' -abar'*ainv];
[v, lambda_eig] = eigs(Adot,Sxxdot,3);
v=v(1:end-1,1)-ainv*v(end,1); %q=1
end

if strcmp(type_quadratic_form,'cov')
S_ab=(a'*b+b'*a); 
abar=sum(a,1)'/(sqrt(NT));
bbar=sum(b,1)'/(sqrt(NT)); 
ainv=pcg(xx,abar,1e-5,1000,Lchol,Lchol');
%binv=pcg(xx,bbar,1e-10,1000,Lchol,Lchol');
Sxxdot=[xx sparse(N+J+K,1) abar+bbar; sparse(1,N+J+K) 1 1; sparse(1,N+J+K) 0 1];
Adot=0.5*[S_ab -S_ab*ainv abar+bbar; bbar' -bbar'*ainv 1; abar' -abar'*ainv 1];
[v, lambda_eig] = eigs(Adot,Sxxdot,3);
v=v(1:end-2,1)-ainv*v(end-1,1); %q=1
end


end

