%function H = hierarchy(A,opts,d)
%-------------------------------------------------------------
%
% Input: A    : adjacency matrix
%        d    : optional positive vector representing the diagonal matrix D = diag(d) 
%        opts : structure with options
%
%-------------------------------------------------------------
% Output: H   : cell containing the hierarchical decomposition of A 
%
% H{j}.cI          : component index of clustering at level j
% H{j}.nc          : number of clusters in level j clustering
% H{j}.A           : diagonally dominant system at level j
% H{j}.invD        : 1./(2*diag(H{j}.A})
% H{j}.dc          : true if diagonal correction occured at level j
% H{j}.islast      : true for the last level of the hierarchy
%
% H{1}.update      : first level directed tree indicating clustering 
%                    [used  for updating hierarchy]
%
% H{last_level}.iterative : true if an iterative method is used 
% 
% H{last_level).chol       : cholesky factorization of H{direct}.A(p,p)
%               chol.ld    : structure containing row major format of ld
%               chol.ldT   : structure containing row major format of ld'
%                   .p     : permutation p
%                   .invp  : inverse permutation of p
%                          
%
%
%c
%----- Options --------------------------------------------           
% opts.update : directed tree indicating an approximate first level
%               clustering. Used to pass the H{1}.update from previous run. 
%
%
% opts.display    : if 1  information is printed (default 0)
%
% opts.direct     : the size which is considered solvable by a direct
%                   method (default 500)
%
%--------------------------------------------------------------

%% 

%/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
%/* 
%/* The CMG solver is distributed under the terms of the GNU General Public   */
%/* Lincense Version 3.0 of the Free Software Foundation.                    */
%/* CMG is also available under other licenses; contact authors for details. */


function H = hierarchy(A,opts,excD) %#ok<FNDEF>



%% input handling

try opts.direct;
catch opts.direct = 500;
end


try opts.display; disp(sprintf('\n'));
catch opts.display =0;
end

try opts.maxit;
catch opts.maxit = 10;
end

%% main hierarchy construction

H{1}.laplacian = 0;

%% reduction of Laplacian to Strictly DD %%
if (nargin < 3)
    n = length(A);
    excD = A(1:(n-1),n);
    A = A(1:(n-1),1:(n-1));
end


TOTAL_TIME = cputime;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% updates temporarily not working
%% try 
%%    updateFlag = length(opts.update);
%% catch
%%    updateFlag = 0;
%%end
updateFlag=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


total_hierarchy_edges = 0;
hierarchy_stagnation =0;
j = 0; lcurrent = length(A); 
while and(lcurrent>opts.direct, j<20)
    
    j = j+1;
    if (j>1)
        A = adjacency_cmg(H{j}.A);       
    end
    
    n = length(A);
    
    % get components cI, row sums DA, and update info for first level
    if (j==1) 
        if (updateFlag == 0)
            [cI, ncomponents, DA, H{1}.update] = steiner_group(A);
        else
            [cI, ncomponents, DA, H{1}.update] = steiner_group(A,opts.update);
        end
    else
        [cI, ncomponents, DA] = steiner_group(A);
    end   
    
    DA = DA + excD;
   


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct level j operators and level j+1 Laplacian and excD

   
    
    if (j==1)
        H{1}.A = laplacian2(A,DA); 
        H{1}.invD = 1./(2*DA);
    else
        H{j}.invD = 1./(2*DA);   
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for strong diagonal dominance
    H{j}.islast = 0;
    H{j}.iterative = 0;
    threshold = (excD>(0.5*DA));
    strong_dominant = find(threshold ==true);  %0.25 is ad-hoc
    not_dominant = find(threshold ==false);
    diag_correction = (length(strong_dominant)>=1);
    H{j}.dc = double(diag_correction);
    if (diag_correction)
       cI(strong_dominant) = uint32(ncomponents+1); 
       H{j}.sd = strong_dominant;

       % check for eliminated clusters and pack index vector
       clear usedI;
       usedI(1:ncomponents,1)  = false;
       usedI(ncomponents+1)    = true; %#ok<AGROW>
       usedI(cI(not_dominant)) = true; %#ok<FNDSB,AGROW>
       ncI = vpack(usedI);
       cI = ncI(cI);
       ncomponents = ncI(length(ncI))-1;
       if (ncomponents == 0)
           H{j}.islast = 1;
           H{j}.iterative = 1;
           hierarchy_stagnation=1; % for good reasons
           break;
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %H{j}.excD = excD;  %for debug

    ambient_ncomponents = double(ncomponents+double(diag_correction));
    Rt = sparse(double(cI),1:n,1,ambient_ncomponents,n);
    R =  Rt';
    H{j+1}.A = Rt*(H{j}.A)*R;
    excD = Rt*excD;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % absorb strong diagonal dominance
    if (diag_correction)
      excD = excD(1:ncomponents);
      excD = excD + abs(H{j+1}.A(1:ncomponents,ambient_ncomponents));
      H{j+1}.A = H{j+1}.A(1:ncomponents,1:ncomponents);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    H{j}.cI = cI -1; % adjust for C-indexing
    H{j}.nc = ambient_ncomponents;
    
    total_hierarchy_edges = total_hierarchy_edges+nnz(H{j+1}.A);
    if or(H{j}.nc >= (length(H{j}.A)-1), total_hierarchy_edges > (4*nnz(H{1}.A)) )
        H{j}.islast =1;
        H{j}.iterative=1;
        hierarchy_stagnation=1; %for potentially bad reasons
        warning('CMG convergence may be slow due to matrix density. Future versions of CMG will eliminate this problem.');
        break;
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lbefore = lcurrent;
    lcurrent = length(H{j+1}.A);     
      
    if (opts.display)
        ff = effective_degrees(adjacency_cmg(H{j+1}.A));
        avgf = sum(ff)/lcurrent;
        min16f = ff(find(ff<(1/8))); %#ok<FNDSB>
        min32f = find(min16f<(1/16));
        min64f = find(min16f<(1/32));
        disp(sprintf('Level %d. Graph vertices: %d ',j+1, lcurrent));
        disp(sprintf('Density : %g',nnz(H{j+1}.A)/length(H{j+1}.A)));
        disp(sprintf('Reduction factor in vertices : %g',lbefore/lcurrent));
        disp(sprintf('Reduction factor in edges    : %g',nnz(H{j}.A)/nnz(H{j+1}.A)));
        disp(sprintf('*****************'));
        %disp(sprintf('Level %d. Graph vertices: %d . Reduction factor:  %g.   Precondition factor average: %g, less than 1/8--1/16--1/32: %g --%g--%g',j+1, lcurrent, lbefore/lcurrent, avgf, length(min16f),length(min32f),length(min64f)  ));
    end

    H{j}.islast = 0; %one level has been added

end %main while loop


if (hierarchy_stagnation == 1)
    last_level = j;
    H{j}.nc = 1;
end

if (hierarchy_stagnation == 0)
last_level = j+1;

if (last_level==1)
    H{last_level}.A = laplacian_cmg(A)+diag(excD);
end

H{last_level}.islast = 1;
H{last_level}.iterative = 0;

prm = amd(H{last_level}.A);
[q,d,ld] = cholGsparse(H{last_level}.A(prm,prm),length(H{last_level}.A));
ld = ld-eye(length(d))+diag(d);
%[ld s prm] = ldlchol(H{last_level}.A);
%ld = (ld - diag(diag(ld))+(diag(diag(ld))).^2)';
H{last_level}.chol.ld = process_sparse_matrix(ld,0,2);
H{last_level}.chol.ldT = process_sparse_matrix(ld',0,2);
H{last_level}.chol.p = uint32(prm-1);
H{last_level}.chol.invp = uint32(inv_permutation(prm)-1);
end % if-hierarchy-stagnation

if (opts.display)
disp(sprintf('Hierarchy constructed in %f seconds',cputime-TOTAL_TIME));
end
%% 
%     Determining recursive calls between levels

if not(H{1}.iterative)
   for k=1:(last_level-1)
       H{k}.repeat = max(floor(nnz(H{k}.A)/nnz(H{k+1}.A)-1),1);
       %H{k}.repeat =1;
       %disp(H{k}.repeat);
   end
end

return