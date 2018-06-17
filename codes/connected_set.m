function [y,id,firmid,id_old,firmid_old,controls,prov_indicator] = connected_set(y,id,firmid,lagfirmid,controls,prov_indicator)
%Finding the largest connected set, using the code from CHK(2013)


%In case we want to carry a province indicator.
if nargin <= 5 
prov_indicator=ones(size(y,1),1);
end 

%Save ids
firmid_old=firmid;
id_old=id;

%RENAME
N=length(y);
sel=~isnan(lagfirmid);

%relabel the firms
[firms,m,n]=unique([firmid;lagfirmid(sel)]);

firmid=n(1:N);
lagfirmid(sel)=n(N+1:end);


%relabel the workers
[ids,m,n]=unique(id);
id=n;

%initial descriptive stats
s=['Original Dataset:'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['# of p-y obs: ' int2str(length(y))];
disp(s);
s=['# of workers: ' int2str(max(id))];
disp(s);
s=['# of firms: ' int2str(max(firmid))];
disp(s);
s=['mean wage: ' num2str(mean(y))];
disp(s)
s=['variance of wage: ' num2str(var(y))];
disp(s)
%FIND CONNECTED SET
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Finding connected set...'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
A=sparse(lagfirmid(sel),firmid(sel),1); %adjacency matrix
%make it square
[m,n]=size(A);
if m>n
    A=[A,zeros(m,m-n)];
end
if m<n
    A=[A;zeros(n-m,n)];
end
A=max(A,A'); %connections are undirected
Adj_matrix=A;

[sindex, sz]=components(A); %get connected sets
idx=find(sz==max(sz)); %find largest set
s=['# of firms: ' int2str(length(A))];
disp(s);
s=['# connected sets:' int2str(length(sz))];
disp(s);
s=['Largest connected set contains ' int2str(max(sz)) ' firms'];
disp(s);

firmlst=find(sindex==idx); %firms in connected set
sel=ismember(firmid,firmlst);

y=y(sel); 
firmid=firmid(sel); 
id=id(sel);
firmid_old=firmid_old(sel);
id_old=id_old(sel);
controls=controls(sel,:);
prov_indicator=prov_indicator(sel,:);


%relabel the firms
[firms,m,n]=unique(firmid);
firmid=n;

%relabel the workers
[ids,m,n]=unique(id);
id=n;
end

