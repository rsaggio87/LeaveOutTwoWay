function elist=index_constr(clustering_var,id,match_id)
NT=length(clustering_var);
count=ones(length(clustering_var),1);
gcs = cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
maxD=max(gcs);
index=(1:NT)';
list_final=[];

for t=1:maxD
    list_base=[clustering_var(gcs==t) index(gcs==t) match_id(gcs==t)];
    LIST_BASE = array2table(list_base,...
    'VariableNames',{'id_cluster','row_count','id'});
        for tt=t:maxD
        list_aux=[clustering_var(gcs==tt) index(gcs==tt)];
        LIST_SEL = array2table(list_aux,...
        'VariableNames',{'id_cluster','column_count'});
        merge=outerjoin(LIST_BASE,LIST_SEL);
        merge = table2array(merge);
        merge=merge(~any(isnan(merge),2),:); %not merged cases are such that a worker does not have future observations. 
        merge=[merge(:,2) merge(:,5) merge(:,3)]; %row, col
        list_final=[list_final;merge];
        end        
end

%Now sort in the appropriate dimension
[~,IX]=sort(list_final(:,1));
elist=list_final(IX,:);
end

