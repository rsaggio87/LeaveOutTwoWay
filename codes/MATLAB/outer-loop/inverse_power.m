%function [V D] = invp(Afun,nAfun,nvec,opts); 
%
% Computes approximations to the largest nvec eigenvectors 
%  of the matrix A, where A*x = Afun(x) for all vectors n
%
% A is assumed to be symmetric
%
% nAfun is the dimension of the argument x of Afun 
%
% opts is structure with the following possible fields
%
% opts.tol: convergence tolerance.
% The iteration stops if two consecutive approximate eigenvectors v,u
% satisfy  1-|v'*u| < tol. Default value tol = 10^(-6)
%
% opts.invpowers is the maximum allowed number of inverse powers iterations 
%   per eigenvector. Default value invpowers = 20;
%
% opts.diplay_info=1 displays convergence information. default is 0.
%
% opts.Vapprox: its columns are used as first approximations of
%  the targeted eigenvectors. 

function [Veig Deig] = invp(Afun,nAfun,nvec,opts); 

%% initializations
warning off;

n = nAfun;
try  tol=opts.tol;
catch tol=10^(-6); end

try invpowers=opts.invpowers; 
catch invpowers=20; end

try display_info=opts.display_info;
catch display_info=0; end

try, 
    Vforce = opts.Vforce;
    Vforceflag = 1;
catch
    Vforceflag = 0; 
end
 
    


%% loop over eigenvectors
for j=1:nvec,
    
       
    % compute initial eigenvector approximation
    if (j==1) 
        try 
            x = opts.Vapprox(:,j); 
        catch 
            x = rand(n,1);
            x = x/norm(x);
        end
    else
        try 
            x = opts.Vapprox(:,j); % use given approximation
        catch
            x = project(x_prior,Veig,1);     % form first approximation from prior approximation
        end
    end

%% loop over inverse powers. computing one eigenvector
     for k=1:invpowers

        if(j==1)
            if (Vforceflag ==1)
                x = project(x,Vforce,1);
            end
            x = x/norm(x);
        else
            if (Vforceflag ==1)
                x = project(x_prior,Vforce,1);
            end
            x = project(x,Veig,1);
        end

        if (k>1)
            if (display_info ==1)
                    disp(sprintf('The last two approximate eigenvectors have angle %e',1-abs(x'*x_prior)));
            end

            if ((1-abs(x'*x_prior)) <tol)
                break;
            end
        end
         
        
        x_prior = x;
        x = Afun(x);
        if (display_info==1)
            y = x/norm(x); 
            rlq = (y'*Afun(y));
            disp(sprintf('The inverse Rayleigh quotient of eigenvector no %d in inverse power %d is %g',j+1,k,rlq));    
        end
        
    end % for
%% end of loop    
    
    if(j==1)
        Veig(:,j) = x/norm(x);
    else
        Veig(:,j) = project(x,Veig,1);
    end
    Deig(j,j) = Veig(:,j)'*Afun(Veig(:,j));

end

      
%warning on