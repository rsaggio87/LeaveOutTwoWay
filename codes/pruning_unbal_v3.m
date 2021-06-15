function [y,firmid,id,id_old,firmid_old,controls,prov_indicator,A] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls,prov_indicator)
% Author: Raffaele Saggio. raffaele.saggio@berkeley.edu

%  This function finds the leave one out largest connected set. 
%  That is, a connected set where if we were to remove the entire working
%  history of a given individual, the resulting graph still remains
%  connected. See appendix B of KSS.
%     
%-INPUTS: a connected bipartite graph where:
% y: outcome variable. Dimensions: N* x 1; N*= # of person-year observations.
% id: worker indicators. Dimensions: N* x 1
% firmid: firm indicators. Dimensions: N* x 1
% id_old: original worker indicators. Dimensions: N* x 1             
% firmid_old: original firm indicators. Dimensions: N* x 1 
% controls: original matrix of controls. Dimensions: N* x K
%     
%  %Output: y, firmid,id, controls in the leave one out sample. firmid, id are all 
%          relabbeled to run estimation in this new sample. firmid_old, 
%          id_old preserve their orginal denomination. 

%          A is the J x J adj matrix associated with leave one out sample
%          where J is the # of Firms in the Leave one out sample. 


if nargin < 7
prov_indicator=ones(size(y,1),1);
end

n_of_bad_workers=1;
while n_of_bad_workers>=1

    %Largest connected set
    A=build_adj(id,firmid); 
    %[sindex, sz]=components(A); %get connected sets
    [sindex,sz] = conncomp(graph(A));
    sz=sz';
    sindex=sindex';
    idx=find(sz==max(sz)); %find largest set
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
    y=y(sel);
    firmid=firmid(sel);
    id=id(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    controls=controls(sel,:);
    prov_indicator=prov_indicator(sel,:);
    
    %Resetting ids one last time.
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;    
    
    
    %Check for leave out connectevity
    NT=size(id,1);
    bad_obs=zeros(NT,1);
    rows=(1:NT)';
    parfor ii=1:NT
    sel=rows~=ii;
    id_sel=id(sel);
    [~,~,id_sel]=unique(id_sel);
    firmid_sel=firmid(sel);
    [~,~,firmid_sel]=unique(firmid_sel);
    A=build_adj(id_sel,firmid_sel);
    [sindex] = conncomp(graph(A)); 
    n_components=unique(sindex);
    bad_obs(ii)=(n_components>1);
    end
    disp('# of observations that when removed would disconnect the graph')
    sum(bad_obs)
    n_of_bad_workers=sum(bad_obs);
    
    %reset
    sel=~bad_obs;
    y=y(sel);
    firmid=firmid(sel);
    id=id(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    controls=controls(sel,:);
    prov_indicator=prov_indicator(sel,:);
    
    %Resetting ids one last time.
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;  
      
end
end

    