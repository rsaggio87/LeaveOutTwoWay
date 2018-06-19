function [y,firmid,id,controls,prov_indicator]=table_1(y,id,firmid,controls,prov_indicator,city)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes as input the person year observations from a particular
%region and then performs the following tasks:

%IMPORTANT: original src data must be sorted by id year (xtset id year)

% 1. Finds and summarizes the largest connected set using routine in CHK (2013).
% 2. Proceeds to find leave out largest connected set, see Appendix B.
% 3. Summarizes leave out largest connected set
% 4. Performs extra pruning for Rovigo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path to matlab bgl


%% STEP 1: FIND LARGEST CONNECTED SET
%Lagfirmid
gcs = [NaN; id(1:end-1)];
gcs = id~=gcs;
lagfirmid=[NaN; firmid(1:end-1)];
lagfirmid(gcs==1)=NaN; %%first obs for each worker

%Largest connected set.
[y,id,firmid,id_old,firmid_old,controls,prov_indicator] = connected_set(y,id,firmid,lagfirmid,controls,prov_indicator);

%Find movers
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
Nmovers=sum(movers);
J=max(firmid);

clear gcs stayer movers T lagfirmid

s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Results on Largest Connected Set '];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(max(Nmovers))];
disp(s);
s=['# of Firms: ' num2str(J)];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);

%% STEP 2: LEAVE ONE OUT CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding the leave one out largest connected set... '];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
tic
[y,firmid,id,~,~,controls,prov_indicator]= pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls,prov_indicator);
disp('Time to find leave one out largest connected set')
toc
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
%Find movers
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
Nmovers=sum(movers);
movers=movers(id);
J=max(firmid);
clear gcs T lagfirmid


%% STEP 3: EXTRA PRUNING FOR ROVIGO ONLY
%We provide one extra pruning for Rovigo to make sure that conditions of
%Theorem 1 are satisfied. In particular, we identify a mover with the
%highest Lindeberg condition and remove the stayers from the two frms this 
%individual moved between. 

%To understand how we identified this particular influential observation,
%comment the following lines and proceed to calculations of the Lindeberg
%condition on the unpruned leave out sample of Rovigo. This will show that 
%mover associated with py observation XXX and XXX is the one with the
%highest Lindeberg condition. 

if strcmp(city,'Rovigo')
firms_bad=[129;151];
sel=ismember(firmid,firms_bad);
sel=and(sel,~movers);
sel=~sel; 
y=y(sel);
id=id(sel);
firmid=firmid(sel);
controls=controls(sel,:);
prov_indicator=prov_indicator(sel,:);     
end

%% STEP 4: SUMMARIZE LEAVE OUT CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['Info on the leave one out connected set:'];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
s=['# of Movers: ' num2str(Nmovers)];
disp(s);
s=['# of Firms: ' num2str(J)];
disp(s);
s=['# of Person Year Observations: ' num2str(size(y,1))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
end

