function invP = inv_permutation(P)

n = length(P);
PM = sparse(P,[1:n],1,n,n);
[invP,P1] = find(PM');



