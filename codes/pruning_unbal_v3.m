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
%find movers
    J=max(firmid);
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
%id and firmids associated with movers only
    id_movers=id(move);
    firmid_mover=firmid(move);

%need to normalize id_movers and keep a dictionary. 
    id_movers_orig=id_movers;
    [~,~,id_movers]=unique(id_movers);
    n_movers=max(id_movers);
    
%unique pairs    
    list=[id_movers firmid_mover];
    list=unique(list,'rows');
    
%Construct Adj. Matrix of Bipartite Graph
    B=sparse(list(:,1),list(:,2),1,n_movers,J);
    G = [ sparse( n_movers, n_movers ), B; B.', sparse(J, J)];
    
%Inspect the resulting Graph: find the workers that constitute an
%articulation point
    [artic_points, CC] = biconnected_components(G);
    bad_workers=artic_points(artic_points<=n_movers);
    
%Now get the right index for these workers
    sel=ismember(id_movers,bad_workers);
    bad_workers=id_movers_orig(sel);
    bad_workers=unique(bad_workers);
    n_of_bad_workers=size(bad_workers,1);
    s=['Number of workers that when dropped would disconnect the graph: ' num2str(n_of_bad_workers)];
    disp(s);
    
%Drop these workers
    sel=~ismember(id,bad_workers);
    y=y(sel);
    firmid=firmid(sel);
    id=id(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    controls=controls(sel,:);
    prov_indicator=prov_indicator(sel,:);
    
%Resetting ids.
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;
    
%Largest connected set once removed bad workers
    A=build_adj(id,firmid); 
    [sindex, sz]=components(A); %get connected sets
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
end
end

    