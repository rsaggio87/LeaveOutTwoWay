function [list_final,nnz_2] = check_clustering(clustering_var)

NT=length(clustering_var);
index=(1:NT)';
[~,~,clustering_var]=unique(clustering_var);
[~,IX]=sort(clustering_var);
clustering_var=clustering_var(IX);
index=index(IX);
count=ones(NT,1);
gcs = cell2mat(accumarray(clustering_var,count,[],@(x){cumsum(x)}));
maxD=max(gcs);
list_final=[];

for t=1:maxD
    list_base=[clustering_var(gcs==t) index(gcs==t)];
    LIST_BASE = array2table(list_base,...
    'VariableNames',{'id_cluster','row_count'});
        for tt=1:maxD
        list_aux=[clustering_var(gcs==tt) index(gcs==tt)];
        LIST_SEL = array2table(list_aux,...
        'VariableNames',{'id_cluster','column_count'});
        merge=outerjoin(LIST_BASE,LIST_SEL);
        merge = table2array(merge);
        merge=merge(~any(isnan(merge),2),:);
        merge=[merge(:,2) merge(:,4) merge(:,2)]; %row, col
        list_final=[list_final;merge];
        end        
end

nnz_2=size(list_final,1);
end

