%function sa = process_sparse_matrix(A,sym,pre)
%
% Input: Sparse matrix A
%      : pre = 2 for double precision, pre =1 for single precision
%      : sym = 1 for saving only symmetric part, 0 otherwise
% Output: Structure containing the row major format of A
%
% 
%

function sa = process_sparse_matrix(A,sym,pre)

    if not(or(sym==1,sym==0))
        error('first input must be 1 or 0');
    end
    
    if not(or(pre==1,pre==2))
        error('second input must be 1 or 2');
    end


n = length(A);
if (sym == 0)
    A=A';
end
if (sym == 1)
    A = tril(A);
end    
    

[i j v] = find(A);

[xx col_j] = unique(j,'first');
col_j = uint32(col_j-1);
col_j(n+1)= uint32(length(v));
row_i = uint32(i-1);

if (pre==1)
  sa.a = single(v);
end
if (pre==2)
  sa.a = v;
end
sa.ia = col_j; 
sa.ja = row_i;
sa.issym = sym;
sa.n = n;

