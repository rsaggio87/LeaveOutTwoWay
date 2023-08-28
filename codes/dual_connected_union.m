function [X_union, X_non_union,NDD] = dual_connected_union(firmid,firmid_orig,peso,union_status,N);


%% reshape things in original NT space
firmid_orig                = repelem(firmid_orig,peso,1); %weight by lenght of the spell
firmid                     = repelem(firmid,peso,1); %weight by lenght of the spell
J                          = max(firmid);

%% connected firms, i.e. firm that have at least one union job and one non-union job
firmid_orig_union_only     = firmid_orig(union_status==1);
firmid_union_only          = firmid(union_status==1);
[~,~,firmid_NORMA]         = unique(firmid_orig_union_only);
firmid_orig_coll	       = accumarray(firmid_NORMA,firmid_orig_union_only,[],@(x)mean(x));
firmid_coll				   = accumarray(firmid_NORMA,firmid_union_only,[],@(x)mean(x));
N_UNION                    = accumarray(firmid_NORMA,1);
elist1					   = [firmid_orig_coll N_UNION  firmid_coll];%% Original Firmid; N of Union Jobs; Id for AKM


firmid_orig_union_only     = firmid_orig(union_status==2);
firmid_union_only          = firmid(union_status==2);
[~,~,firmid_NORMA]         = unique(firmid_orig_union_only);
firmid_orig_coll	       = accumarray(firmid_NORMA,firmid_orig_union_only,[],@(x)mean(x));
firmid_coll				   = accumarray(firmid_NORMA,firmid_union_only,[],@(x)mean(x));
N_UNION                    = accumarray(firmid_NORMA,1);
elist2					   = [firmid_orig_coll N_UNION firmid_coll];%% Original Firmid; N of Non-Union Jobs;  Id for AKM

LIST_BASE 				   = array2table(elist1,...
    								'VariableNames',{'firmid_original','N_union','ID_AKM_UNION'});
LIST_SEL 				   = array2table(elist2,...
        							'VariableNames',{'firmid_original','N_NON_union','ID_AKM_NON_UNION'});
merge					   = outerjoin(LIST_BASE,LIST_SEL);
merge 					   = table2array(merge);
merge					   = merge(~any(isnan(merge),2),:); %not merged cases are firms with either 100% of union workers or 0% 
peso                       = merge(:,2)+merge(:,5);
dual_list                  = [merge(:,1) peso merge(:,3) merge(:,6)]; %original firmid; counts of union + non union workers; firmid if looking at union x firm effect; firmid if looking at non-union x firm effect;
NDD                        = size(dual_list,1);

ids                        = sparse((1:NDD)',dual_list(:,3)',1,NDD,J);
S                          = speye(J-1);
S                          = [S;sparse(-zeros(1,J-1))];
ids                        = [ids*S];
X_union                    = [sparse(NDD,N) ids]; %number and ids of columns as in the original AKM
X_union                    = repelem(X_union,peso,1); %% so variance component is weighted by number of person-year obs for a given firm 

ids                        = sparse((1:NDD)',dual_list(:,4)',1,NDD,J);
S                          = speye(J-1);
S                          = [S;sparse(-zeros(1,J-1))];
ids                        = [ids*S];
X_non_union                = [sparse(NDD,N) ids]; %number and ids of columns as in the original AKM
X_non_union                = repelem(X_non_union,peso,1); %% so variance component is weighted by number of person-year obs for a given firm 

end