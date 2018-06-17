function [elist,index_unique,row,column,Tweight] = input_lambda_P_FE_fast(firmid_delta,firmid_delta_f,id_movers,Tweight)
%this function keeps track of the non zero elements along the (block) diagonal 
%of the hat matrix  X(X'*X)^(-1)X'. We need to keep track of this when 
%computing the (matrix) of statistical leverages xi'(X'*X)^(-1)xi for each
%individual i.


sel=firmid_delta~=firmid_delta_f;
index=(1:size(id_movers,1))';
D_movers=sparse(index(sel),id_movers(sel)',1,size(id_movers,1),max(id_movers));
D_movers=D_movers*D_movers';
[index1, index2]=find(triu(D_movers));
clear D_movers
elist=[firmid_delta(index1) firmid_delta_f(index1) firmid_delta(index2) firmid_delta_f(index2)];
row=index1;
column=index2;
[elist, ~, index_unique]=unique(elist,'rows');
Tweight=Tweight(index1);
end

