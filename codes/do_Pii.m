function Lambda_P = do_Pii(X,clustering_var)
n=size(X,1);

%If no clustering, Lambda_P is just diagonal matrix.
if isempty(clustering_var)
	clustering_var = (1:n)';
end

%Set matrices for parallel environment. 
xx=X'*X;
xx_c = parallel.pool.Constant(xx);
X_c =  parallel.pool.Constant(X);
Lchol=ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
Chol_c =  parallel.pool.Constant(Lchol);

%Options for the solver. 
numIterations = 300; %iteration for the pcg solver
tol=1e-6; %tol for pcg


%Return the structure of the indexes associated with the clustering variable
elist = check_clustering(clustering_var);
M=size(elist,1);

%Set elist in parallel environment.
elist_1=parallel.pool.Constant(elist(:,1));
elist_2=parallel.pool.Constant(elist(:,2));
Pii=zeros(M,1);

parfor i=1:M
    [xtilde, flag]= pcg(xx_c.Value,X_c.Value(elist_2.Value(i),:)',tol,numIterations,Chol_c.Value,Chol_c.Value');
    Pii(i)=X_c.Value(elist_2.Value(i),:)*xtilde;
end        

Lambda_P=sparse(elist(:,1),elist(:,2),Pii,n,n);

end

