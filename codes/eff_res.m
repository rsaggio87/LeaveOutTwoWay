function [Lambda_P, Lambda_B_fe, Lambda_B_cov, Lambda_B_pe] = eff_res(X,xx,Lchol,N,J,K,elist,clustering_level,movers,T,type_algorithm,id,firmid,epsilon,tolProb)
%This function calculates, using parallel coding, (Pii,Bii) for a general 
%two way fixed effects model. The code is likely to be slow on large
%datasets.

%Read
do_cov=1;
do_pe=1;
if nargout<2
    error('More output should be specified')
end
if nargout<=2
    do_cov=0;
    do_pe=0;
end
if nargout<=3
    do_pe=0;
end

if nargin<=13
    tolProb=0.5;
    epsilon=0.01;
end
if nargin<=14
    tolProb=0.5;
end

%Dimensions
NT=size(X,1);
M=size(elist,1);

%PreCreate
Pii=zeros(M,1);
Bii_fe=zeros(M,1);
Bii_cov=zeros(M,1);
Bii_pe=zeros(M,1);

%Options for solver
numIterations = 300; %iteration for the pcg solver
tol=1e-6; %tol for pcg

%Objects for parfor 
    xx_c = parallel.pool.Constant(xx);
    X_c =  parallel.pool.Constant(X);
    Xright=parallel.pool.Constant(X(elist(:,2),:));
    Xleft=parallel.pool.Constant(X(elist(:,1),:));
    if K > 0
        Lchol_c=parallel.pool.Constant(Lchol);
    end

%Special case of Laplacian design matrix    
if K == 0
    [~,~,match_id]=unique([id firmid],'rows');
    Nmatches=max(match_id);
    match_id_movers=match_id(movers);
    firmid_movers=firmid(movers);
    id_movers=id(movers);
    [~,sel,~]=unique(match_id_movers,'rows');
    match_id_movers=match_id_movers(sel);
    firmid_movers=firmid_movers(sel);
    id_movers=id_movers(sel);
    maxT=max(T(~movers));
    count=ones(NT,1);
    gcs = cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
    sel_stayers=(gcs==1).*(~movers);
    sel_stayers=sel_stayers>0;
    stayers_matches_sel=match_id(sel_stayers);
    Tinv=1./T;
    elist_JLL=[id_movers N+firmid_movers id_movers N+firmid_movers];
    M=size(elist_JLL,1);
    Pii_movers=zeros(M,1);
    Bii_fe_movers=zeros(M,1);
    Bii_cov_movers=zeros(M,1);
    Bii_pe_movers=zeros(M,1);
end
    
%Loop
if K>0
    if ~strcmp(clustering_level,'obs')
        parfor i=1:M     
            [xtilde_right, flag]= pcg(xx_c.Value,Xright.Value(i,:)',tol,numIterations,Lchol_c,Lchol_c.Value');
            [xtilde_left, flag]= pcg(xx_c.Value,Xleft.Value(i,:)',tol,numIterations,Lchol_c.Value,Lchol_c.Value');

            %Statistical Leverage
            Pii(i)=Xleft.Value(i,:)*xtilde_right;

            %Bii for Variance of Firm Effects
            aux_right=xtilde_right(N+1:N+J-1,:);
            aux_left=xtilde_left(N+1:N+J-1,:);
            COV=cov(X_c.Value(:,N+1:N+J-1)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
            Bii_fe(i)=COV(1,2)*(NT-1);

            %Bii for Variance of Person Effects
            if do_pe == 1
                aux_right=xtilde_right(1:N);
                aux_left=xtilde_left(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,1:N)*aux_right);
                Bii_pe(i)=COV(1,2)*(NT-1);
            end

            %Bii for Covariance of Person, Firm Effects
            if do_cov == 1
                aux_right=xtilde_right(N+1:N+J-1);
                aux_left=xtilde_left(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
                Bii_cov(i)=COV(1,2)*(NT-1);
            end
        end
    end

    if strcmp(clustering_level,'obs')
        parfor i=1:M     
        [xtilde, flag]= pcg(xx_c.Value,Xright.Value(i,:)',tol,numIterations,Lchol_c.Value,Lchol_c.Value');

        %Statistical Leverage
        Pii(i)=Xright.Value(i,:)*xtilde;

        %Bii for Variance of Firm Effects
        aux_right=xtilde(N+1:N+J-1,:);
        aux_left=xtilde(N+1:N+J-1,:);
        COV=cov(X_c.Value(:,N+1:N+J-1)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
        Bii_fe(i)=COV(1,2)*(NT-1);


        %Bii for Variance of Person Effects
            if do_pe==1
                aux_right=xtilde(1:N);
                aux_left=xtilde(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,1:N)*aux_right);
                Bii_pe(i)=COV(1,2)*(NT-1);
            end

         %Bii for Covariance of Person, Firm Effects
            if do_cov==1
                aux_right=xtilde(N+1:N+J-1);
                aux_left=xtilde(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
                Bii_cov(i)=COV(1,2)*(NT-1);
            end
        end
    end
end

%Loop
if K==0
    
    if  strcmp(type_algorithm,'exact')
        Xright = sparse((1:M)',elist_JLL(:,1),1,M,N+J);
        Xright =  Xright+sparse((1:M)',elist_JLL(:,2),-1,M,N+J);
        Xright=parallel.pool.Constant(Xright);
        
        parfor i=1:M
        [xtilde, flag]= pcg(xx_c.Value,Xright.Value(i,:)',tol,numIterations,Lchol);

        %Statistical Leverage
        Pii_movers(i)=Xright.Value(i,:)*xtilde;

        %Bii for Variance of Firm Effects
        aux_right=xtilde(N+1:N+J-1,:);
        aux_left=xtilde(N+1:N+J-1,:);
        COV=cov(X_c.Value(:,N+1:N+J-1)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
        Bii_fe_movers(i)=COV(1,2)*(NT-1);


        %Bii for Variance of Person Effects
            if do_pe==1
                aux_right=xtilde(1:N);
                aux_left=xtilde(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,1:N)*aux_right);
                Bii_pe_movers(i)=COV(1,2)*(NT-1);
            end

         %Bii for Covariance of Person, Firm Effects
            if do_cov==1
                aux_right=xtilde(N+1:N+J-1);
                aux_left=xtilde(1:N);
                COV=cov(X_c.Value(:,1:N)*aux_left,X_c.Value(:,N+1:N+J-1)*aux_right);
                Bii_cov_movers(i)=COV(1,2)*(NT-1);
            end 
        end
        
    end
    
   if  strcmp(type_algorithm,'JLL') 
       
       %Number of random draws to implement Random Projection Algorithm.
       disp('# of Simulated Projections for JLL:')
       
       %scale = ceil(log(N+J)/epsilon^2); Actual bound provided by S-S
       %(2011). 
       
       %In my experience, the formula above is too conservarive. One can obtain very
       %good approximation in very reasonable computation time using the
       %following alternative found here: http://www.cs.cmu.edu/~jkoutis/EffectiveResistances/EffectiveResistances.m
       
       scale = ceil(log2(NT))/epsilon
       
       %Auxiliary components
       X_fe=parallel.pool.Constant([sparse(NT,N) X(:,N+1:end)]); 
       X_pe=parallel.pool.Constant([X(:,1:N) sparse(NT,J)]);
       elist_1=parallel.pool.Constant(elist_JLL(:,1));
       elist_2=parallel.pool.Constant(elist_JLL(:,2));
       elist_3=parallel.pool.Constant(elist_JLL(:,3));
       elist_4=parallel.pool.Constant(elist_JLL(:,4));        
       parfor i=1:scale
                
                %To avoid redundant warnings from each worker.
                Z=0;
                ZB=0;
                ZB_pe=0;
                
                %Random Projection Matrix (Rademacher)
                ons = (rand(1,NT) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                
                %Lambda_P
                [Z, flag]= pcg(xx_c.Value,(ons*(X_c.Value))',tol,numIterations,Lchol);
                
                %Lambda_B_fe
                ons=ons-mean(ons);
                [ZB, flag]= pcg(xx_c.Value,(ons*X_fe.Value)',tol,numIterations,Lchol);
                
                %Lambda_B_cov
                if do_cov == 1
                    [ZB_pe, flag]= pcg(xx_c.Value,(ons*X_pe.Value)',tol,numIterations,Lchol);
                end
                
                %Collect results of the given draw
                Pii_movers=Pii_movers+(((Z(elist_1.Value)-Z(elist_2.Value))).*(((Z(elist_3.Value)-Z(elist_4.Value)))));
                Bii_fe_movers=Bii_fe_movers+(((ZB(elist_1.Value)-ZB(elist_2.Value))).*(((ZB(elist_3.Value)-ZB(elist_4.Value)))));
                
                if do_cov == 1
                    Bii_cov_movers=Bii_cov_movers+(((ZB(elist_1.Value)-ZB(elist_2.Value))).*(((ZB_pe(elist_3.Value)-ZB_pe(elist_4.Value)))));
                end
                
                if do_pe == 1
                    Bii_pe_movers=Bii_pe_movers+(((ZB_pe(elist_1.Value)-ZB_pe(elist_2.Value))).*(((ZB_pe(elist_3.Value)-ZB_pe(elist_4.Value)))));
                end
       end
                
   end
  
   
   %Assign step
   Pii_movers=sparse(match_id_movers,1,Pii_movers,Nmatches,1);
   Pii_stayers=sparse(stayers_matches_sel,1,Tinv(sel_stayers),Nmatches,1);
   Pii=Pii_movers+Pii_stayers;
   
   Bii_fe=sparse(match_id_movers,1,Bii_fe_movers,Nmatches,1);
   
   if do_cov == 1
        Bii_cov=sparse(match_id_movers,1,Bii_cov_movers,Nmatches,1);
   end
   
   if do_pe  == 1
        Bii_pe=sparse(match_id_movers,1,Bii_pe_movers,Nmatches,1);
        stayers=~movers;
         for t=2:maxT %T=1 have Pii=1 so need to be dropped.
            sel=(gcs==1).*stayers.*(T==t);
            sel=sel>0;
            index_sel=find(sel);
            match_sel_aux=match_id(sel);
            first=index_sel(1);
            Xuse=X(first,:);
            [xtilde, flag]= pcg(xx,Xuse',tol,numIterations,Lchol);
            aux_right=xtilde(1:N);
            aux_left=xtilde(1:N);
            COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
            Bii_pe_stayers=COV(1,2)*(NT-1);
            Bii_pe_stayers=sparse(match_sel_aux,1,Bii_pe_stayers,Nmatches,1);
            Bii_pe=Bii_pe+Bii_pe_stayers;
         end    
   end
   
end

%Create the matrices.
rows=elist(:,1);
column=elist(:,2);
index_cluster=elist(:,3);

if K == 0 
    Pii=Pii(index_cluster);
    Bii_fe=Bii_fe(index_cluster);
    if do_cov==1
        Bii_cov=Bii_cov(index_cluster);
    end
    if do_pe==1
        Bii_pe=Bii_pe(index_cluster);
    end
end

%Lambda P
Lambda_P=sparse(rows,column,Pii,NT,NT);
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.

%Lambda B var(fe)
Lambda_B_fe=sparse(rows,column,Bii_fe,NT,NT);
Lambda_B_fe=Lambda_B_fe+triu(Lambda_B_fe,1)'; %make it symmetric.

%Lambda B cov(fe,pe)
if do_cov==1
    Lambda_B_cov=sparse(rows,column,Bii_cov,NT,NT);
    Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)';
end

%Lambda B, var(pe)
if do_pe==1
    Lambda_B_pe=sparse(rows,column,Bii_pe,NT,NT);
    Lambda_B_pe=Lambda_B_pe+triu(Lambda_B_pe,1)'; %make it symmetric.
end

end

