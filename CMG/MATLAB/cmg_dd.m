% function pfun = cmg_dd(A,opts)
%
% A     : adjacency matrix (positive off-diagonals, zero diagonal)
% opts  : structure with options
%
% opts.d : optional non-negative vector 
%
% pfun : function implementing a preconditioner for 
%        the matrix laplacian(A)+opts.d
%
% 
%
% For an arbitrary b-side in the null space of A, 
%    the system Ax = b can be solved by 
%    x = pcg(A, b, tol,iter, pfun);        % >> help pcg for documentation



% update: directed forest indicating previously computed first level clustering
% update:  directed forest indicating new first level clustering



function [pfun H] = cmg_dd(A,opts) %#ok<STOUT,FNDEF>

default_lap=1;

opts.display = 1;
precision = 1;

try
    
    if not(length(A)==length(opts.d))
      error('Dimensions of opts.d and A must agree');
    end   
    
    if (min(opts.d)<0)
       error('opts.d must be a non-negative vector');
    end


    lap = 0;
    s = sum(opts.d);
    sA = sum(A); 
    smax = max(max(sA),s);
    smin = min(min(sA),s);
    
    if ((smin/smax)<1e-07)
        precision =2; 
    end    
    
    if ((smin/smax)<1e-015)
         precision =2;
         disp('Warning: The system is numerically ill-conditioned.')
    end  
    

catch %#ok<CTCH>
    lap = 1;
    sA = sum(A); 
    if ((min(sA)/max(sA))<1e-07)
        precision =2; 
    end

    if ((max(sA)/min(sA))<1e-015)
         precision =2;
         disp('Warning: The system is numerically ill-conditioned.'); 
    end    
    
end

precision=2;
if (precision == 1)
if (lap<1)
   H = hierarchy(A,opts,opts.d);
   %stats = hierarchy_stat_analysis(H);
   H{1}.laplacian = lap;
   sH = process_hierarchy(H,1,1);
   pfun = @(b) double(mx_s_preconditioner(sH,single(b)));
end

if (lap==1)
   H = hierarchy(A,opts); 
   H{1}.laplacian = lap;
   sH = process_hierarchy(H,1,1);
   pfun = @(b) double(mx_s_preconditioner(sH,single(b)));
end
end

if (precision == 2)
if (lap<1)
   H = hierarchy(A,opts,opts.d);
   %stats = hierarchy_stat_analysis(H);
   H{1}.laplacian = lap;
   sH = process_hierarchy(H,1,2);
   pfun = @(b) mx_d_preconditioner(sH,b);
end

if (lap==1)
   H = hierarchy(A,opts); 
   H{1}.laplacian = lap;
   sH = process_hierarchy(H,1,2);
   pfun = @(b) mx_d_preconditioner(sH,b);
end
end




