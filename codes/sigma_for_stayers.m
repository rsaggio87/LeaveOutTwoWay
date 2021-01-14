function sigma_stayers = sigma_for_stayers(y,id,firmid,peso,b);
%Computes heteroskedatic robust variances for stayers 

%Get back to person-year space
id=repelem(id,peso,1);
firmid=repelem(firmid,peso,1);

%Compute Pii for stayers (1/T_i)
T=accumarray(id,1);
Pii=1./T;
Mii=1-Pii(id);

%Compute OLS residual
NT=size(y,1);
D=sparse(1:NT,id',1);
F=sparse(1:NT,firmid',1);
J=size(F,2);
S= speye(J-1);
S= [S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
X= [D,F*S];
eta=y-X*b;
eta_h=eta./Mii; %Leave one out residual
sigma_stayers=(y-mean(y)).*eta_h;

%Collapse back to match
[~,~,match_id] 		= unique([id firmid],'rows','stable');
sigma_stayers       = accumarray(match_id,sigma_stayers,[],@(x)mean(x));       
end

