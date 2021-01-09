function [Pii, Mii, correction_JLA, Bii_fe, Bii_cov, Bii_pe] = leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,scale)
%This function calculates, using parallel coding, (Pii,Bii) for a general 
%two way fixed effects model and the unbiased estimate of
%Var(e_gt) where e_gt is the error term of the AKM model

%Read
do_cov=1;
do_pe=1;
tolProb=0.5;

if nargout<4
    error('More output should be specified')
end

if nargout==4
    do_cov=0;
    do_pe=0;
end

if nargout==5
    do_pe=0;
end


%Dimensions
NT=size(X,1);

%PreCreate
Pii     = sparse(NT,1);
Bii_fe  = sparse(NT,1);
Bii_cov = sparse(NT,1);
Bii_pe  = sparse(NT,1);
Mii	 	= sparse(NT,1);
Pii_sq  = sparse(NT,1);
Mii_sq  = sparse(NT,1);
Pii_Mii = sparse(NT,1);

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
        disp('Running Exact Algorithm...')
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
  correction_JLA=1;
  Mii=1-Pii;
 end

%Random Projection
 if strcmp(type_algorithm,'JLA')

	%Clean memory
	Lchol    = [];
	xx		 = [];
	X		 = [];

	%disp('# of Simulated Projections used for JLA Algorithm:')
	%scale

        disp('Running JLA Algorithm...')
        parfor i=1:scale     
		%To avoid redundant warnings from each worker.
		Z=0;
		Z_fe=0;
		Z_pe=0;
		Z_cov=0;
		
		%Random Projection Matrix (Rademacher)
		ons 		= (rand(1,NT) > tolProb);
		ons 		= ons - not(ons);
		%ons 		= ons./sqrt(scale);
		
		%Get me the row
		[Z, flag]	= pcg(xx_c.Value,(ons*(X_c.Value))',tol,numIterations,Lchol_c.Value);
		Z		 	= X_c.Value*Z;
		
		%Collect (construct augmented estimator that combines Mii,Pii)
		aux		 	= Z.^2/scale;
		Pii	 	 	= Pii+aux;
		aux			= Z.^4/scale;
		Pii_sq		= Pii_sq + aux;
		
		aux			= ((ons'-Z).^2)/scale;
		Mii			= Mii+aux;
		aux			= ((ons'-Z).^4)/scale;
		Mii_sq		= Mii_sq + aux;
		
		Pii_Mii		= Pii_Mii + ((Z.^2).*((ons'-Z).^2))/scale;
		
		aux			= [];
		
		%Random Projection Matrix (Rademacher)
		ons 		= (rand(1,size(X_fe.Value,1)) > tolProb);
		ons 		= ons - not(ons);
		ons 		= ons./sqrt(scale);
		ons			= ons-mean(ons);

        %Bii for Variance of Firm Effects
        [Z flag]	= pcg(xx_c.Value,(ons*(X_fe.Value))',tol,numIterations,Lchol_c.Value);
        Z_fe		= X_c.Value*Z;				
	    Bii_fe		= Bii_fe+(Z_fe.*Z_fe);	
        
	    %Bii for Variance of Person Effects
        [Z flag]	= pcg(xx_c.Value,(ons*(X_pe.Value))',tol,numIterations,Lchol_c.Value);
        Z_pe		= X_c.Value*Z;				
	    Bii_pe		= Bii_pe+(Z_pe.*Z_pe);
	    
	    %Bii for CoVariance of Person,Firm Effects
        Z_cov		= (Z_pe.*Z_fe);			
	    Bii_cov		= Bii_cov+Z_cov;
        end
        
%A potential issue when estimating Pii with the JLA is that it is possible 
%that such estimate>1. To avoid such issues, which are typically present in
%moderate size datasets with relatively low mobility, we use an estimator
%of 1-Pii that is guaranteed to be in [0,1] and that has smaller variance

%We also account for the nonlinear bias introduced by the JLA procedure
%using the following step. See appendix of Di Addario, Kline, Saggio and
%Soelvsten (2020) for  details. 

Pii                 = Pii./(Pii+Mii); 	
Mii                 = 1-Pii;
Vi                  = (1/scale)*((Mii.^2).*Pii_sq+(Pii.^2).*Mii_sq-2*Mii.*Pii.*Pii_Mii);
Bi                  = (1/scale)*(Mii.*Pii_sq-Pii.*Mii_sq+2*(Mii-Pii).*Pii_Mii);
correction_JLA 		= (1-Vi./(Mii.^2)+Bi./Mii);        
end

 

end

