function [ydelta, Fdelta, Tweight,id_delta,firmid_delta,firmid_delta_f,gcs_delta,list_final]= stacked_Fdelta(y,id,firmid,gcs)
%this function takes the original dataset (expressed in the typical
%person-year format) and creates a new dataset where all the (unique)
%within person time differneces are stacked together.


T=accumarray(id,1);
maxT=max(T);
T=T(id);
list_cell=cell(maxT);

parfor t=1:maxT
    list_base=[firmid(gcs==t) id(gcs==t) gcs(gcs==t) y(gcs==t) T(gcs==t)];
    LIST_BASE = array2table(list_base,...
    'VariableNames',{'Firm_t','id','t','y_t','Ti'});
        list_tt=[];
        for tt=t+1:maxT
        list_aux=[firmid(gcs==tt) id(gcs==tt) gcs(gcs==tt) y(gcs==tt)];
        LIST_SEL = array2table(list_aux,...
        'VariableNames',{'Firm_tprime','id','tprime','y_tprime'});
        merge=outerjoin(LIST_BASE,LIST_SEL);
        merge = table2array(merge);
        merge=merge(~any(isnan(merge),2),:); %not merged cases are such that a worker does not have future observations. 
        merge=[merge(:,1) merge(:,6) merge(:,2) merge(:,3) merge(:,8) merge(:,4) merge(:,9) merge(:,5)]; %Firm_t, Firm_t' ID t t' y_t y_t' Ti
        list_tt=[list_tt;merge];
        end
list_cell{t}=list_tt;
end
list_cell=cell2mat(list_cell);
[~,IX]=sort(list_cell(:,3));
list_final=list_cell(IX,:);

%Now sparsify the list to obtain Fdelta.
Tweight=list_final(:,8);
m=size(list_final,1);
Fdelta  = sparse([1:m 1:m],[list_final(:,1) list_final(:,2)],[-1./sqrt(Tweight) 1./sqrt(Tweight)],m,max(firmid)); %GLS version.
ydelta=list_final(:,7)-list_final(:,6);
ydelta = ydelta./sqrt(Tweight); %already weighted.
count=ones(length(ydelta),1);
gcs_delta = cell2mat(accumarray(list_final(:,3),count,[],@(x){cumsum(x)})); %%this provides a within person count of the differences.
firmid_delta=list_final(:,1);
firmid_delta_f=list_final(:,2);
id_delta=list_final(:,3);
list_final=[list_final gcs_delta];



%ydelta: PD of the outcome variable (Weighted) 
%Fdelta: PD of the matrix of Firm assignments (weighted).
%id_movers: report the id of the worker in the PD space.
%firmid_delta: firmid of firm in period t
%firmid_delta_f: firmid of firm in period t'
%gcs_delta: cumulative count of the pairwise differences.

end

