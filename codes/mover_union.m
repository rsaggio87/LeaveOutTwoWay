function [X_union, X_non_union,NDD] = mover_union(firmid,id,id_orig,union_status);

J                          = max(firmid);
N                          = max(id);


%% connected firms, i.e. firm that have at least one union job and one non-union job
id_orig_union_only     = id_orig(union_status==1);
id_union_only          = id(union_status==1);
[~,~,id_NORMA]         = unique(id_orig_union_only);
id_orig_coll	       = accumarray(id_NORMA,id_orig_union_only,[],@(x)mean(x));
id_coll				   = accumarray(id_NORMA,id_union_only,[],@(x)mean(x));
N_UNION                = accumarray(id_NORMA,1);
elist1				   = [id_orig_coll N_UNION  id_coll];%% Original ID; N of YEARS in Union; Id for AKM


id_orig_union_only     = id_orig(union_status==0);
id_union_only          = id(union_status==0);
[~,~,id_NORMA]         = unique(id_orig_union_only);
id_orig_coll	       = accumarray(id_NORMA,id_orig_union_only,[],@(x)mean(x));
id_coll				   = accumarray(id_NORMA,id_union_only,[],@(x)mean(x));
N_UNION                = accumarray(id_NORMA,1);
elist2				   = [id_orig_coll N_UNION  id_coll];%% Original ID; N of YEARS in Union; Id for AKM
LIST_BASE 				   = array2table(elist1,...
    								'VariableNames',{'id_original','N_union','ID_AKM_UNION'});
LIST_SEL 				   = array2table(elist2,...
        							'VariableNames',{'id_original','N_NON_union','ID_AKM_NON_UNION'});
merge					   = outerjoin(LIST_BASE,LIST_SEL);
merge 					   = table2array(merge);
merge					   = merge(~any(isnan(merge),2),:); %not merged cases are non-switchers
peso                       = merge(:,2)+merge(:,5);
dual_list                  = [merge(:,1) peso merge(:,3) merge(:,6)]; %original id; T; worker id if looking at union x worker id effect; id if looking at non-union x worker id effect; T_non_union
disp('Number of Workers switching union status:')
NDD                        = size(dual_list,1)

ids                        = sparse((1:NDD)',dual_list(:,3)',1,NDD,N);
F                          = sparse(NDD,J-1);
X_union                    = [ids F]; %number and ids of columns as in the original AKM
X_union                    = repelem(X_union,peso,1); %% so variance component is weighted by number of person-year obs for a given worker id 

ids                        = sparse((1:NDD)',dual_list(:,4)',1,NDD,N);
F                          = sparse(NDD,J-1);
X_non_union                = [ids F]; %number and ids of columns as in the original AKM
X_non_union                = repelem(X_non_union,peso,1); %% so variance component is weighted by number of person-year obs for a given worker id 
end