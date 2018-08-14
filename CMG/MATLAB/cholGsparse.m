% [Q,D,L] = cholGsparse(A,k)
%
% Q: Shcur complement of A with respect to the elimination of the
% first k vertices
% D: k x k diagonal
% L: lower triangular such that: A = L[D;0|0;Q]L^T

function [Q,D,L] = cholG(A,r)

n = length(A);
L = speye(n);

for I=1:r
    vi=A(I+1:n,I); di=A(I,I); D(I)=di;
    L(I+1:n,I) = vi/di;
    Z = vi*vi'/di;
    [Iz Jz V] = find(Z); 
    Z = sparse(Iz+I,Jz+I,V,n,n);
    A = A -Z;
end

Q = A(r+1:n,r+1:n);










    

