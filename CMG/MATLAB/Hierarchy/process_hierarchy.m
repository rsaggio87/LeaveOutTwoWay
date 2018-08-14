% function dH = process_hierarchy(H,sym,pre)
%
% Input  H: Output of function hierarchy
%      sym: 1 to exploit symmetry, 0 otherwise
%      pre: 2 for double precision, 1 for single precision
%
% Output dH: Processed hierarchy

function H = process_hierarchy(H,sym,pre)

nH = max(size(H));

if (nH>1)
if (H{nH-1}.islast == 1)
    nH = nH -1; 
end
end

for j=1:(nH)
    
    if not(or(sym==1,sym==0))
        error('first input must be 1 or 0');
    end
    
    if not(or(pre==1,pre==2))
        error('second input must be 1 or 2');
    end
        
    if and(pre == 1,or(j<nH,H{j}.iterative==1))
       H{j}.invD = single(H{j}.invD); 
    end
    
 
    if  or(j<nH,H{j}.iterative==1)
    if  (pre == 1)
        H{j}.lws1 = single(zeros(length(H{j}.A)+1,1));
        H{j}.lws2 = single(zeros(length(H{j}.A)+1,1));
        H{j}.sws1 = single(zeros(H{j}.nc+1,1));
        H{j}.sws2 = single(zeros(H{j}.nc+1,1));
        H{j}.sws3 = single(zeros(H{j}.nc+1,1));        
    else
        H{j}.lws1 = zeros(length(H{j}.A)+1,1);
        H{j}.lws2 = zeros(length(H{j}.A)+1,1);
        H{j}.sws1 = zeros(H{j}.nc+1,1);
        H{j}.sws2 = zeros(H{j}.nc+1,1);
        H{j}.sws3 = zeros(H{j}.nc+1,1);        
    end
    end
   
    H{j}.A = process_sparse_matrix(H{j}.A,sym,pre);
    
    



    
    % add update

    
end

