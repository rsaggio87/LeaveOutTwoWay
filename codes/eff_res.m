function [Lambda_P, Lambda_B_fe, Lambda_B_cov, Lambda_B_pe] = eff_res(X,xx,Lchol,N,J,K,elist,clustering_level,movers,T)
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

%slice
    Xright=X(elist(:,2),:);
    Xleft=X(elist(:,1),:);

if K == 0 && strcmp(clustering_level,'obs')  
    Xright=X(movers,:);
    Nmovers=sum(movers);
    maxT=max(T(~movers));
    movers_index=find(movers);
    Pii_movers=zeros(Nmovers,1);
    Bii_fe_movers=zeros(Nmovers,1);
    Bii_cov_movers=zeros(Nmovers,1);
    Bii_pe_movers=zeros(Nmovers,1);
    Tinv=1./T;
end

%Loop
if K>0
    if ~strcmp(clustering_level,'obs')
        parfor i=1:M     
            [xtilde_right, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol,Lchol');
            [xtilde_left, flag]= pcg(xx,Xleft(i,:)',tol,numIterations,Lchol,Lchol');

            %Statistical Leverage
            Pii(i)=Xleft(i,:)*xtilde_right;

            %Bii for Variance of Firm Effects
            aux_right=xtilde_right(N+1:N+J-1,:);
            aux_left=xtilde_left(N+1:N+J-1,:);
            COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
            Bii_fe(i)=COV(1,2)*(NT-1);

            %Bii for Variance of Person Effects
            if do_pe == 1
                aux_right=xtilde_right(1:N);
                aux_left=xtilde_left(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
                Bii_pe(i)=COV(1,2)*(NT-1);
            end

            %Bii for Covariance of Person, Firm Effects
            if do_cov == 1
                aux_right=xtilde_right(N+1:N+J-1);
                aux_left=xtilde_left(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
                Bii_cov(i)=COV(1,2)*(NT-1);
            end
        end
    end

    if strcmp(clustering_level,'obs')
        parfor i=1:M     
        [xtilde, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol,Lchol');

        %Statistical Leverage
        Pii(i)=Xright(i,:)*xtilde;

        %Bii for Variance of Firm Effects
        aux_right=xtilde(N+1:N+J-1,:);
        aux_left=xtilde(N+1:N+J-1,:);
        COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
        Bii_fe(i)=COV(1,2)*(NT-1);


        %Bii for Variance of Person Effects
            if do_pe==1
                aux_right=xtilde(1:N);
                aux_left=xtilde(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
                Bii_pe(i)=COV(1,2)*(NT-1);
            end

         %Bii for Covariance of Person, Firm Effects
            if do_cov==1
                aux_right=xtilde(N+1:N+J-1);
                aux_left=xtilde(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
                Bii_cov(i)=COV(1,2)*(NT-1);
            end
        end
    end
end

%Loop
if K==0
    if ~strcmp(clustering_level,'obs')
        parfor i=1:M     
            [xtilde_right, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol);
            [xtilde_left, flag]= pcg(xx,Xleft(i,:)',tol,numIterations,Lchol);

            %Statistical Leverage
            Pii(i)=Xleft(i,:)*xtilde_right;

            %Bii for Variance of Firm Effects
            aux_right=xtilde_right(N+1:N+J-1,:);
            aux_left=xtilde_left(N+1:N+J-1,:);
            COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
            Bii_fe(i)=COV(1,2)*(NT-1);

            %Bii for Variance of Person Effects
            if do_pe == 1
                aux_right=xtilde_right(1:N);
                aux_left=xtilde_left(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
                Bii_pe(i)=COV(1,2)*(NT-1);
            end

            %Bii for Covariance of Person, Firm Effects
            if do_cov == 1
                aux_right=xtilde_right(N+1:N+J-1);
                aux_left=xtilde_left(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
                Bii_cov(i)=COV(1,2)*(NT-1);
            end
        end
    end

    if strcmp(clustering_level,'obs')
        parfor i=1:Nmovers     
        [xtilde, flag]= pcg(xx,Xright(i,:)',tol,numIterations,Lchol);

        %Statistical Leverage
        Pii_movers(i)=Xright(i,:)*xtilde;

        %Bii for Variance of Firm Effects
        aux_right=xtilde(N+1:N+J-1,:);
        aux_left=xtilde(N+1:N+J-1,:);
        COV=cov(X(:,N+1:N+J-1)*aux_left,X(:,N+1:N+J-1)*aux_right);
        Bii_fe_movers(i)=COV(1,2)*(NT-1);


        %Bii for Variance of Person Effects
            if do_pe==1
                aux_right=xtilde(1:N);
                aux_left=xtilde(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
                Bii_pe_movers(i)=COV(1,2)*(NT-1);
            end

         %Bii for Covariance of Person, Firm Effects
            if do_cov==1
                aux_right=xtilde(N+1:N+J-1);
                aux_left=xtilde(1:N);
                COV=cov(X(:,1:N)*aux_left,X(:,N+1:N+J-1)*aux_right);
                Bii_cov_movers(i)=COV(1,2)*(NT-1);
            end
        end
        
        %Now use tricks to assign to movers vs. stayers.
        Bii_fe=sparse(movers_index,1,Bii_fe_movers,NT,1);
        Pii=sparse(movers_index,1,Pii_movers,NT,1);
        Pii(~movers)=Tinv(~movers);
        
        if do_cov==1
            Bii_cov=sparse(movers_index,1,Bii_cov_movers,NT,1);
        end
        
        if do_pe==1
            Bii_pe=sparse(movers_index,1,Bii_pe_movers,NT,1);
            stayers=~movers;
            for t=2:maxT %T=1 have Pii=1 so need to be dropped.
            sel=stayers.*(T==t);
            index_sel=find(sel);
            Xuse=X(index_sel,:);  
            Xuse=Xuse(1,:);
            [xtilde, flag]= pcg(xx,Xuse',tol,numIterations,Lchol);
            aux_right=xtilde(1:N);
            aux_left=xtilde(1:N);
            COV=cov(X(:,1:N)*aux_left,X(:,1:N)*aux_right);
            Bii_pe_stayers=COV(1,2)*(NT-1);
            Bii_pe(index_sel)=Bii_pe_stayers;
            end
        end
    end

    
end

%Create the matrices.
rows=elist(:,1);
column=elist(:,2);
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

