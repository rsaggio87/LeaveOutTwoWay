function [test_statistic, linear_combination, SE_linear_combination_KSS, SE_linear_combination_RES,  SE_linear_combination_NAI, stat, pvalue] = lincom_KSS(y,X,Z,Transform,clustering_var,Lambda_P,labels,restrict,NSIM)
%----------------------------------------------------
%%%Description
%----------------------------------------------------
% This function performs inference on linear combinations of 
% regression coefficients suitable for settings with many regressors  
% and potentially heteroscedastic errors. 

% For theoretical justification and further details, see Proposition 1 and 
% Remark 9 of Kline, Saggio, Solvsten (2018) - KSS henceforth.

% We also provide an option that allows inference when the errors are dependent 
% within "clusters." This option requires that the model remains estimable after 
% any given cluster is left out. 

% The linear combinations of interest are specified as v'beta where 
% beta=(X'X)X'y is an OLS coefficient vector of dimension K x 1
% and v is a user specified MATRIX of dimension K x r. 

% To construct v, we assume the linear combinations of interest can be 
% written as projection coefficients, so that: v'=(Z'Z)^(-1)Z'Transform:

%  			- "Z" is a user specified matrix of dimension p x r. 
% 			- "Transform" is a user specified matrix of dimension p x K 

% "Transform" determines which elements of beta=(X'X)X'y are to be projected 
% onto the basis defined by Z. See example below for intuition.

% The function also provides a test that all of the specified 
% linear combinations are jointly equal to zero. This test is based
% upon the quadratic form beta'*A*beta with the matrix A specified as
%
%                  A=(1/r)*v*(v'(X'X)^(-1)v)^(-1)v', 

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
% worker age, firm size, firm age. 

% Let Transform=[sparse(n,N) F] where F is a n x J matrix that collects
% firm indicators. 

% Note that Transform*beta gives back the firm effects. 

% Therefore here the r x 1 vector v'*beta returns the linear regression
% coefficients of a projection of firm effects (person-year weighted) onto
% a constant (automatically added by the function) and 
% the variables contained in the matrix Z. 

% "lincom_KSS" outputs standard errors, t-stat and a joint significant test
% of these r linear combinations.
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
%cluster_var: This determines whether user wants to conduct cluster robust inference.
			 
%			  If cluster_var = [], inference will be robust only to many regressors asymptotics and heteroskedasticity.
			  
%			  If cluster_var ~= [], inference will be robust to many regressors asymptotics, heteroskedasticity and clustering
%			  where the level of clustering is specified by the variable "cluster_var" itself. 
%			  For instance, if cluster_var=worker_id, then the function will output standard errors, t-statistics, etc that 
%			  account for serial correlation in the error term within each worker. 
%----------------------------------------------------
%Lambda_P: a n x n (sparse) matrix that contains the statistical leverages
%          of the high dimensional model, i.e Lambda_P(i,j)=X(i,:)*(X'*X)^(-1)X(j,:)'.

%		   Note that, depending on the value "cluster_var", Lambda_P shall either be
%          a diagonal matrix (when cluster_var =[]) or a block diagonal matrix 
%		   (when cluster_var ~=[]). 
%          
%          This is a non-mandatory input. If the user does not have calculated
%          the (cluster)-leave-out leverages in a prior step, "lincom_KSS" 
%          will compute them in accordance with the mandatory input "cluster_var".
%          (e.g. if cluster_var= worker_id, "lincom_KSS" will calculate the leverages 
%          by leaving out the entire set of observations belonging to a worker) 
		 
%		   When beta refers to a coefficient vector coming from a two-way fixed effects (AKM) model, 
%	 	   we highly recommend to compute Lambda_P in a prior step by calling the function "leave_out_COMPLETE" 
%	 	   which is able to compute the statistical leverages in a very efficient way. Then
%	 	   invoke "lincom_KSS" making sure to specify "Lambda_P" as an input, see examples
% 		   in "example_testing".

%  		   Keep in mind that if Lambda_P is specified as an input, its structure must be coherent 
% 		   with the level of clustering specified by the input "cluster_var". For instance if "cluster_var"=worker_id
%          but "Lambda_P" is a diagonal matrix (which would imply only heteroskedastic robust inference) the code
%		   will produce an error.
%		
%----------------------------------------------------
%labels: non-mandatory input. This is used to label the r columns of Z. It
%        should be provided as a cell.
%----------------------------------------------------
%restrict: non-mandatory input. This should be a matrix that is read only
%          in the case in which the user does not want to conduct a standard 
%          joint test for the joint significance of the columns of Z but 
%          instead wants to calculate a different joint hypothesis 
%          with new vector v defined this time as v' = restrict * v'. See
%          example 4 in "example_testing" for more intuition on how to
%          specify this input.
%----------------------------------------------------
%Nsim: # of simulations to derive the p-value of the joint statistic.
%      Default is 10K.
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
%							  When the clustering option is turned on, the code automatically
%							  augments SE_linear_combination_KSS to account for serial correlation 
%							  of the error term within the dimension specified by the user
%							  in the input cluster_var.
%----------------------------------------------------
%  SE_linear_combination_RES: This is a vector of dimension r x 1 and returns the 
%						      standard error of each linear contrast embedded in v'*beta 
%     						  replacing the term \hat{\sigma}_{i}^{2} in KSS remark 9 with
%							  \tilde{\sigma}_{i}^2=(y-X*beta).^2.

%							  When the clustering option is turned on, the code automatically
%							  augments SE_linear_combination_RES to account for serial correlation 
%							  of the error term within the dimension specified by the user
%							  in the input cluster_var.
%----------------------------------------------------
%----------------------------------------------------
% SE_linear_combination_NAI:  This is a vector of dimension r x 1 and returns the standard
%							  error that one would obtain without accounting for the noise
%							  contained in beta.
%							  
%							  For instance:
%							  in the example writte above, "SE_linear_combination_NAI" will 
%							  give back the standard error that one would obtain when the user
%                             in a first step obtains the firm effects (firm_effects) from an AKM
% 							  model and then in second step types in Stata: "reg firm_effects gender age firm_size firm_age, robust"   
%     					
%                             When the clustering option is turned on, the code automatically
%							  augments "SE_linear_combination_NAI" to account for serial correlation 
%							  of the error term within the dimension specified by the user
%							  in the input cluster_var, i.e. it returns:
%								
%							  "reg firm_effects gender age firm_size firm_age, cluster(cluster_var)" 	
%----------------------------------------------------
%----------------------------------------------------
%  stat: This is a scalar that reports the KSS adjusted quadratic form beta'*A*beta.

%  In particular, stat reports the statistic associated with the null hypothesis that 
%  all the columns in Z jointly predict Transform*beta.

%  			For instance, in the example above when Z = [gender age firm_age firm_size]: 

%  			stat reports the statistic that tests the joint significance of
%  			all the variables contained in Z in predicting the firm effects. 
%  			p_value reports the associated p-value obtained by means of simulations.
%  
%  See example 4 in "example_testing" to undestand how one can specify an
%  additional matrix "restrict" to conduct more flexible joint tests.
%----------------------------------------------------
%  p_value: p-value associated stat. Value obtained via simulations.
%----------------------------------------------------

%% READ INPUTS
got_Pii=0;
got_labels=0;
got_restrict=0;
if nargin < 9 
    NSIM=10000;
end
if nargin==6 && ~isempty(Lambda_P)
    got_Pii=1;
end
if nargin==7 && ~isempty(Lambda_P) && ~isempty(labels)
    got_Pii=1;
    got_labels=1;
end
if nargin==7 && isempty(Lambda_P) && ~isempty(labels)
    got_labels=1;
end
if nargin==7 && ~isempty(Lambda_P) && isempty(labels)
    got_Pii=1;
end
if nargin==8 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end
if nargin==8 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_labels=1;
    got_restrict=1;
end
if nargin==8 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_restrict=1;
end
if nargin==8 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==8 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==8 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end
if nargin==9 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end

if nargin==9 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict) 
    got_labels=1;
    got_restrict=1;
end
if nargin==9 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict) 
    got_Pii=1;
    got_restrict=1;
end
if nargin==9 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==9 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==9 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end

if nargin==10 && ~isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict)
    got_Pii=1;
    got_labels=1;
    got_restrict=1;
end

if nargin==10 && isempty(Lambda_P) && ~isempty(labels) && ~isempty(restrict) 
    got_labels=1;
    got_restrict=1;
end
if nargin==10 && ~isempty(Lambda_P) && isempty(labels) && ~isempty(restrict) 
    got_Pii=1;
    got_restrict=1;
end
if nargin==10 && ~isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_Pii=1;
    got_labels=1;
end 
if nargin==10 && isempty(Lambda_P) && ~isempty(labels) && isempty(restrict)
    got_labels=1;
end
if nargin==10 && ~isempty(Lambda_P) && isempty(labels) && isempty(restrict)
    got_Pii=1;
end

%% SET DIMENSIONS
n=size(X,1);
K=size(X,2);
%% Add Constant
Z=[ones(size(Z,1),1) Z];
%% PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
xx=X'*X;
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
xy=X'*y;
[beta flag]=pcg(xx,xy,1e-10,1000,Lchol,Lchol');
eta=y-X*beta;

%% PART 1B: VERIFY LEAVE OUT COMPUTATION
if got_Pii == 0 
    Lambda_P=do_Pii(X,clustering_var);
end

if got_Pii == 1 && ~isempty(clustering_var)
    nnz_1=nnz(Lambda_P);
    [~,nnz_2]=check_clustering(clustering_var);
    if nnz_1 == nnz_2
    	s=['******************************************'];
    	disp(s);
		check=['The structure of the specified Lambda_P is consistent with the level of clustering required by the user.'];
		disp(check)
		s=['******************************************'];
    	disp(s);
	end
	if nnz_1 ~= nnz_2
		error('The user wants cluster robust inference but the Lambda_P provided by the user is not consistent with the level of clustering asked by the user. Try to omit input Lambda_P when running lincom_KSS' )
	end
end
I_Lambda_P=(speye(n,n)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
msg = lastwarn ; 

if ~isempty(strfind(msg, 'singular'))
error('Cluster-Leave Out Residual cannot be computed. Possibly, this occurred because the level of clustering specified does not allow to compute the leave out OLS coefficient. Example: user sets cluster_var=worker_id but beta includes worker effects--> code cannot identify a worker effect after leaving out all the observations related to that worker')
end

%% PART 2: SET UP MATRIX FOR SANDWICH FORMULA
[rows,columns]=find(Lambda_P); %read the non zero elements of Lambda_P 
aux= 0.5*(y(rows).*eta_h(columns)+y(columns).*eta_h(rows));
sigma_i=sparse(rows,columns,aux,n,n);
aux= 0.5*(eta(rows).*eta(columns)+eta(columns).*eta(rows));
sigma_i_res=sparse(rows,columns,aux,n,n);
r=size(Z,2);
wy=Transform*beta;
zz=Z'*Z;
numerator=Z\wy;
chet=wy-Z*numerator;
aux= 0.5*(chet(rows).*chet(columns)+chet(columns).*chet(rows));
sigma_i_chet=sparse(rows,columns,aux,n,n);


%% PART 3: COMPUTE
denominator=zeros(r,1);
denominator_RES=zeros(r,1);
for q=1:r
    v=sparse(q,1,1,r,1);
    v=zz\v;
    v=Z*v;
    v=Transform'*v;
    [right flag]=pcg(xx,v,1e-5,1000,Lchol,Lchol');
    left=right';    
    denominator(q)=left*(X'*sigma_i*X)*right;
    denominator_RES(q)=left*(X'*sigma_i_res*X)*right;
end    
test_statistic=numerator./(sqrt(denominator));
zz_inv=zz^(-1);
SE_linear_combination_NAI=zz_inv*(Z'*sigma_i_chet*Z)*zz_inv;

%% PART 4: REPORT
	s=['******************************************'];
    disp(s);
    disp(s);
    disp('INFERENCE ON LINEAR COMBINATIONS')
    s=['******************************************'];
    disp(s);  
    disp(s); 
if got_labels == 0 
    for q=2:r
    if q <= r    
    s=['Linear Combination - Column Number ' num2str(q-1) ' of Z: ' num2str(numerator(q))];
    disp(s)
    s=['Standard Error of the Linear Combination - Column Number ' num2str(q-1) ' of Z: ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat - for Column Number ' num2str(q-1) ' of Z: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s);
    end
    end
end

if got_labels == 1 
    for q=2:r
    tell_me=labels{q-1}; 
    s=['Linear Combination associated with '  tell_me ':  ' num2str(numerator(q))];
    disp(s)
    s=['SE Linear Combination associated with '  tell_me ':  ' num2str(sqrt(denominator(q)))];
    disp(s)
    s=['T-stat: ' num2str(test_statistic(q))];
    disp(s)
    s=['******************************************'];
    disp(s)
    end   
end

%% PART 5: Joint-test. Quadratic form beta'*A*beta
if nargout > 5
    if got_restrict == 0 
        restrict=sparse((1:r-1)',2:r,1,r-1,r); %omitting the constant from the test.
    end
    v=restrict*(zz\(Z'*Transform));
    v=v';
    v=sparse(v);
    r=size(v,2); 
    
%Auxiliary
   aux=xx\v;
   opt_weight=v'*aux;
   opt_weight=opt_weight^(-1);
   opt_weight=(1/r)*(opt_weight+opt_weight')/2; %force it to be symmetric for eigs to work properly.
   
%Eigenvalues, eigenvectors, and relevant components.
    FunAtimesX = @(x) v*opt_weight*(v'*x);
    opts.issym = 1; 
    opts.tol = 1e-10;
    [V, lambda]=eigs(FunAtimesX,K,xx,r,'largestabs',opts);
    lambda=diag(lambda);
    W=X*V;
    V_b=W'*sigma_i*W;    
    
%Now focus on obtaining matrix Lambda_B with the A test associated with a joint hypothesis testing.
    Bii=full(opt_weight)^(1/2)*aux'; %r x k
    Bii=Bii*X'; %r x n
    Bii=Bii'; %n x r;
    Bii = 0.5*(Bii(rows,:).*Bii(columns,:) + Bii(columns,:).*Bii(rows,:)); 
    Bii = sum(Bii,2);
    Lambda_B=sparse(rows,columns,Bii,n,n);
    
%Leave Out Joint-Statistic
    stat=(v'*beta)'*opt_weight*(v'*beta)-y'*Lambda_B*eta_h;
    
%Now simulate critical values under the null.
    mu=zeros(r,1);    
    sigma = V_b;
    b_sim = mvnrnd(mu,sigma,NSIM);
    theta_star_sim=sum(lambda'.*(b_sim.^2 - diag(V_b)'),2); %This only works appropriately with Matlab2017
    pvalue=mean(theta_star_sim>stat);  
    
%Report
    s=['******************************************'];
    disp(s);
    disp(s);
    disp('JOINT TEST')
    s=['******************************************'];
    disp(s);
    disp(s);
    s=['Joint-Test Statistic:   ' num2str(stat)];
    disp(s)
    s=['p-value:    ' num2str(pvalue)];
    disp(s)
    s=['******************************************'];
    disp(s)
    disp(s)
    
    
end
test_statistic=test_statistic(2:end);
linear_combination=numerator(2:end);
SE_linear_combination_KSS=sqrt(denominator(2:end));
SE_linear_combination_RES=sqrt(denominator_RES(2:end));
SE_linear_combination_NAI=diag(SE_linear_combination_NAI);
SE_linear_combination_NAI=sqrt(SE_linear_combination_NAI(2:end));
end


