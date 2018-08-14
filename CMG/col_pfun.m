function M = col_pfun( A, pfun )

[n, m] = size(A);
for k=1:m
    M(:,k) = pfun(A(:,k)); 
    M(:,k) = M(:,k)-sum(M(:,k))/n;
end

end

