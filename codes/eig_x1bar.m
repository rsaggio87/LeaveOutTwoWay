function [lambda,x1bar] = eig_x1bar(X,Q,lambda_eig,trace)
lambda=(lambda_eig.^2)/trace;
q=Q(:,1); %the corresponding eigenvector.
%calculate x1bar
x1bar=X*q;
norm=(sum(x1bar.^2))^(0.5);
x1bar=x1bar/norm;
end

