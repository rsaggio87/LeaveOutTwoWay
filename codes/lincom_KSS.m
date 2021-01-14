function [test_statistic, linear_combination, SE_linear_combination_KSS] = lincom_KSS(y,X,Z,Transform,sigma_i,labels)
%----------------------------------------------------
%%%Description
%----------------------------------------------------
% This function performs inference on linear combinations of 
% regression coefficients suitable for settings with many regressors  
% and potentially heteroscedastic errors. 

% For theoretical justification and further details, see Proposition 1 and 
% Remark 9 of Kline, Saggio, Solvsten (2018) - KSS henceforth.


% The linear combinations of interest are specified as v'beta where 
% beta=(X'X)X'y is an OLS coefficient vector of dimension K x 1
% and v is a user specified MATRIX of dimension K x r. 

% To construct v, we assume the linear combinations of interest can be 
% written as projection coefficients, so that: v'=(Z'Z)^(-1)Z'Transform:

%  			- "Z" is a user specified matrix of dimension p x r. 
% 			- "Transform" is a user specified matrix of dimension p x K 

% "Transform" determines which elements of beta=(X'X)X'y are to be projected 
% onto the basis defined by Z. See example below for intuition.

% NOTE: this function is appropriate for cases where r is "small" meaning
% that it does not grow with the sample size asymptotically (See Proposition 1 of KSS).
% Moreover, this function will only work with versions of Matlab2017A or higher.

%----------------------------------------------------
%%%Example
%----------------------------------------------------
% Let y be a n x 1 vector of wages.

% Let X be a n x K "AKM" matrix collecting 
% worker and firm indicators. There are  a total of N worker effects and J  
% firm effects so that K=N+J.

% Let Z be a n x r matrix of "extra" variables such as worker gender,
% worker age, firm size, firm age, region indicators. 

% Let Transform=[sparse(n,N) F] where F is a n x J matrix that collects
% firm indicators. 

% Note that Transform*beta gives back the firm effects. 

% Therefore here the r x 1 vector v'*beta returns the linear regression
% coefficients of a projection of firm effects (person-year weighted) onto
% a constant (automatically added by the function) and 
% the variables contained in the matrix Z. 

% "lincom_KSS" outputs standard errors and t-stat for each of the columns
% present in Z
%----------------------------------------------------
%%%Inputs
%----------------------------------------------------
%y: outcome vector. Dimension n x 1.
%----------------------------------------------------
%X: regressor matrix. Dimension n x k.
%----------------------------------------------------
%Z: matrix defining the linear contrast. Dimension p x r.
%	Note: A constant will always be added to Z by the code.
%----------------------------------------------------
%Transform: matrix that selects / weights elements of the coefficient vector beta. Dimension p x k.
%----------------------------------------------------
%sigma_i: KSS (leave-out) estimate of the heteroskedastic variances present
%in the original model
%		
%----------------------------------------------------
%labels: non-mandatory input. This is used to label the r columns of Z. It
%        should be provided as a cell.

%----------------------------------------------------
%----------------------------------------------------
%%%Outputs
%----------------------------------------------------
%----------------------------------------------------
%  test_statistic: This is a vector of dimension r x 1. 
%  
%  Each element of this vector reports the t-stat associated with one of 
%  the r linear combination embedded in v.

%  			For instance, in the example above, let Z = [gender age firm_age firm_size]. 

%  			Then: 
%  			test_statistic(1) reports the t-statistic on the gender coefficient from 
%  			a regression of firm effects onto a constant, a gender dummy and firm
%  			size.

%  			test_statistic(2) reports the t-statistic on the age coefficient from 
%  			a regression of firm effects onto a constant, a gender dummy and firm
%  			size.

%  			test_statistic(3) reports the t-statistic on the firm_age coefficient from 
%  			a regression of firm effects onto a constant, a gender dummy and firm
%  			size.

%  			test_statistic(4) reports the t-statistic on the firm_size coefficient from 
%  			a regression of firm effects onto a constant, a gender dummy and firm
%  			size.
%----------------------------------------------------
%  linear_combination: This is a vector of dimension r x 1. This returns v'*beta.
%----------------------------------------------------
%----------------------------------------------------
%  SE_linear_combination_KSS: This is a vector of dimension r x 1 and returns the 
%						      standard error of each linear contrast embedded in v'*beta 
%     						  as derived in KSS remark 9.
%							
%							  Note that if the function lincom_KSS is called 
%                             within the function "leave_out_KSS" then the
%                             reported standard errors are going to be robust
%                             to serial correlation within match.
%% READ INPUTS
got_labels=0;
if nargin==6 && ~isempty(labels) 
    got_labels=1;
end

%% SET DIMENSIONS
n                   = size(X,1);
r                   = size(Z,2)+1;
%% Add Constant
Z                   =[ones(size(Z,1),1) Z];
%% COMPUTE
xx                  = X'*X;
Lchol               = ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy                  = X'*y;
beta                = pcg(xx,xy,1e-10,1000,Lchol,Lchol');
wy                  = Transform*beta;
zz                  = Z'*Z;
numerator           = Z\wy;
sigma_i             = sparse((1:n)',1:n,sigma_i,n,n);
sigma_i_naive       = sparse((1:n)',1:n,(y-X*beta).^2,n,n);
denominator         = zeros(r,1);
denominator_naive   = zeros(r,1);
for q=1:r
    v=sparse(q,1,1,r,1);
    v=zz\v;
    v=Z*v;
    v=Transform'*v;
    [right flag]=pcg(xx,v,1e-5,1000,Lchol,Lchol');
    left=right';    
    denominator(q)=left*(X'*sigma_i*X)*right;
    denominator_naive(q)=left*(X'*sigma_i_naive*X)*right;
end    
test_statistic=numerator./(sqrt(denominator));

%% PART 3: REPORT
	s=['******************************************'];
    disp(s);
    disp(s);
    disp('RESULTS ON LINCOM')
    s=['******************************************'];
    disp(s);  
    disp(s); 
if got_labels == 0 
    for q=2:r
    if q <= r    
    s=['Coefficient on - Column Number ' num2str(q-1) ' of Z: ' num2str(numerator(q))];
    disp(s)
    s=['Robust "White" Standard Error: ' num2str(sqrt(denominator_naive(q)))];
    disp(s)
    s=['KSS Standard Error: ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s);
    end
    end
end

if got_labels == 1 
    for q=2:r
    tell_me=labels{q-1}; 
    s=['Coefficient on '  tell_me ':  ' num2str(numerator(q))];
    disp(s)
    s=['Robust "White" Standard Error: ' num2str(sqrt(denominator_naive(q)))];
    disp(s)
    s=['KSS Standard error:  ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s)
    end   
end
test_statistic=test_statistic(2:end);
linear_combination=numerator(2:end);
SE_linear_combination_KSS=sqrt(denominator(2:end));
end


