function A=build_adj(id,firmid)
gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    stayer=(firmid==lagfirmid);
    stayer(gcs==1)=1;
    stayer=accumarray(id,stayer);
    T=accumarray(id,1);
    stayer=T==stayer;
    movers=stayer~=1;
    movers=movers(id);
    id_movers=id(movers);
    id_movers_orig=id_movers;
    [ids,m,id_movers]=unique(id_movers);
    firmid_movers=firmid(movers);
    gcs_mover=gcs(movers);
    lagfirmid=[NaN; firmid_movers(1:end-1)];
    lagfirmid(gcs_mover==1)=NaN; %%first obs for each worker
    list=[lagfirmid firmid_movers id_movers_orig]; %%all the moves observed from one period to the next this time only for movers
    sel=~isnan(list(:,1));
    list=list(sel,:);
    sel=list(:,1)~=list(:,2);
    list=list(sel,:);
    [e,jj,kk]=unique(list(:,1:2),'rows');
    weight=histc(kk,1:numel(jj)); %frequencies
    J=max(firmid);
    A = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[weight;weight],J,J); %adjacency matrix
end

