function [theta,theta_KSS,sigma_i,Lambda_P, Lambda_B] = kss_computation(y,D,F,A_1,A_2,scale)

numIterations = 10000; %iteration for the pcg solver
tol=1e-10; %tol for pcg
NT=size(y,1);
X=[D F];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section 1: PI ESTIMATES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xy                                  = X'*y;
[beta,coeff_DF,xx,Lchol]            = areg_KSS(xy,D,F);
eta                                 = y-X*beta;
right                               = A_2*beta;
left                                = A_1*beta;
theta                               = cov(left,right);
theta                               = theta(1,2);
beta=[];
xy=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section 2: GETTING THINGS READY TO RUN JLA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PreCreate
Pii	 	= zeros(NT,1);
Mii	 	= zeros(NT,1);
Bii	 	= zeros(NT,1);
Pii_sq  = zeros(NT,1);
Mii_sq  = zeros(NT,1);
Pii_Mii = zeros(NT,1);
same = isequal(A_1,A_2);

%Options for solver of JLA
numIterations = 1000; %iteration for the pcg solver
tol=1e-6; %tol for pcg
tolProb=0.5;

%Objects for parfor 
xx_c 	=  parallel.pool.Constant(xx);
X_c 	=  parallel.pool.Constant(X);
Lchol_c =  parallel.pool.Constant(Lchol);
A_1 	=  parallel.pool.Constant(A_1);
A_2 	=  parallel.pool.Constant(A_2);
D       =  parallel.pool.Constant(D);
F       =  parallel.pool.Constant(F);
coeff_DF= parallel.pool.Constant(coeff_DF);


%Clean memory
Lchol   = [];
xx		= [];
X		= [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section 3: RUN JLA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
		%ons 		= ons./sqrt(scale);
		
		%Get me the row
        Z           = areg_KSS((ons*(X_c.Value))',D.Value,F.Value,coeff_DF.Value,xx_c.Value,Lchol_c.Value);
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
		ons 		= (rand(1,size(A_1.Value,1)) > tolProb);
		ons 		= ons - not(ons);
		ons 		= ons./sqrt(scale);
		ons			= ons-mean(ons); %we always assume A=A_1'*Q*A_2 where Q is a symmetric demeaning matrix

        %Bii
        [Z_2 flag]	= areg_KSS((ons*(A_2.Value))',D.Value,F.Value,coeff_DF.Value,xx_c.Value,Lchol_c.Value);
        Z_2			= X_c.Value*Z_2;
        if same == 1	
        	Z_1		= Z_2;
        end
        if same == 0	
        	[Z_1 flag]	= areg_KSS((ons*(A_1.Value))',D.Value,F.Value,coeff_DF.Value,xx_c.Value,Lchol_c.Value);
        	Z_1			= X_c.Value*Z_1;	
        end
        
        Z			= (Z_1.*Z_2);			
	    Bii			= Bii+Z;
	    
        end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section 4: COMBINE JLA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	
	Pii = Pii./(Pii+Mii); 	
	Mii = 1-Pii;
	Vi  = (1/scale)*((Mii.^2).*Pii_sq+(Pii.^2).*Mii_sq-2*Mii.*Pii.*Pii_Mii);
	Bi  = (1/scale)*(Mii.*Pii_sq-Pii.*Mii_sq+2*(Mii-Pii).*Pii_Mii);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section 5: COMPUTE KSS ESTIMATOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta_h				= eta./Mii; %Leave one out residual
sigma_i				= (y-mean(y)).*eta_h; %KSS estimate of individual variance, we de-meaned the outcome here.
dof					= size(left,1)-1;
correction_JLA 		= (1-Vi./(Mii.^2)+Bi./Mii);
sigma_i				= sigma_i.*correction_JLA; %to adjust for non-linear bias
theta_KSS			= theta-(1/dof)*sum(Bii.*sigma_i);

%Lambda P
Lambda_P			= sparse(1:NT,(1:NT)',Pii,NT,NT);

%Lambda B
Lambda_B			= sparse(1:NT,(1:NT)',Bii,NT,NT);


end

