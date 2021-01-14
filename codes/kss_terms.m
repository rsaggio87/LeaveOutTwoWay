function [Pii, correction_JLA,Bii_fe, Bii_cov, Bii_pe] = leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,scale)
%This function calculates, using parallel coding, (Pii,Bii) for a general 
%two way fixed effects model and the unbiased estimate of
%Var(e_gt) where e_gt is the error term of the AKM model

%Read
do_cov=1;
do_pe=1;
tolProb=0.5;

if nargout<3
    error('More output should be specified')
end

if nargout<=4
    do_cov=0;
    do_pe=0;
end

if nargout<=4
    do_pe=0;
end

%Dimensions
NT=size(X,1);

%PreCreate
Pii=zeros(NT,1);
Bii_fe=zeros(NT,1);
Bii_cov=zeros(NT,1);
Bii_pe=zeros(NT,1);

%Options for solver
numIterations = 300; %iteration for the pcg solver
tol=1e-6; %tol for pcg

%Objects for parfor 
xx_c 	=  parallel.pool.Constant(xx);
X_c 	=  parallel.pool.Constant(X);
Lchol_c =  parallel.pool.Constant(@() Lchol);
X_fe 	=  parallel.pool.Constant(X_fe);
X_pe 	=  parallel.pool.Constant(X_pe);
    
%Exact calculations 
 if strcmp(type_algorithm,'exact')
        parfor i=1:NT     
        [xtilde, flag]= pcg(xx_c.Value,X_c.Value(i,:)',tol,numIterations,Lchol_c.Value);

        %Statistical Leverage
        Pii(i)=X_c.Value(i,:)*xtilde;

        %Bii for Variance of Firm Effects
       	COV=cov(X_fe.Value*xtilde,X_fe.Value*xtilde);
        Bii_fe(i)=COV(1,2)*(size(X_fe.Value,1)-1);


        %Bii for Variance of Person Effects
        if do_pe==1
        		COV=cov(X_pe.Value*xtilde,X_pe.Value*xtilde);
        		Bii_pe(i)=COV(1,2)*(size(X_pe.Value,1)-1);
        end

         %Bii for Covariance of Person, Firm Effects
        if do_cov==1
        		COV=cov(X_pe.Value*xtilde,X_fe.Value*xtilde);
        		Bii_cov(i)=COV(1,2)*(size(X_pe.Value,1)-1);
        end
        end
 end

%Random Projection
 if strcmp(type_algorithm,'JLL')

	%Clean memory
	Lchol    = [];
	xx		 = [];
	X		 = [];

	disp('# of Simulated Projections for JLL:')
	scale

 
        parfor i=1:scale     
		%To avoid redundant warnings from each worker.
		Z=0;
		Z_fe=0;
		Z_pe=0;
		Z_cov=0;
		
		%Random Projection Matrix (Rademacher)
		ons 		= (rand(1,NT) > tolProb);
		ons 		= ons - not(ons);
		ons 		= ons./sqrt(scale);
		
		%Get me the row
		[Z, flag]	= pcg(xx_c.Value,(ons*(X_c.Value))',tol,numIterations,Lchol_c.Value);
		
		%Collect
		Z		 	= X_c.Value*Z;
		Z		 	= Z.^2;
		Pii	 	 	= Pii+Z;
		
		%Random Projection Matrix (Rademacher)
		ons 		= (rand(1,size(X_fe.Value,1)) > tolProb);
		ons 		= ons - not(ons);
		ons 		= ons./sqrt(scale);
		ons			= ons-mean(ons);

        %Bii for Variance of Firm Effects
        [Z flag]	= pcg(xx_c.Value,(ons*(X_fe.Value))',tol,numIterations,Lchol_c.Value);
        Z			= X_c.Value*Z;	
        Z_fe		= (Z.*Z);			
	    Bii_fe		= Bii_fe+Z_fe;
	    
	    %Bii for Variance of Person Effects
        [Z flag]	= pcg(xx_c.Value,(ons*(X_pe.Value))',tol,numIterations,Lchol_c.Value);
        Z			= X_c.Value*Z;	
        Z_pe		= (Z.*Z);			
	    Bii_pe		= Bii_pe+Z_pe;
	    
	    %Bii for CoVariance of Person,Firm Effects
        Z_cov		= (Z_pe.*Z_fe);			
	    Bii_cov		= Bii_cov+Z_cov;
 		end
 end		

%Obtain the estimate of heteroskedastic variance of each error term

eta_h				= eta./Mii; %Leave one out residual
sigma_i				= (y-mean(y)).*eta_h; %KSS estimate of individual variance, we de-meaned the outcome here.
dof					= size(left,1)-1;
correction_JLA 		= (1-Vi./(Mii.^2)+Bi./Mii);
sigma_i				= sigma_i.*correction_JLA; %to adjust for non-linear bias
 
 

end

