function [y, D, clusterID, year, X] = create_DGP(wins_value,N_clusters,N_X,T,seed_s)
%% start by creating a balanced panel
clusterID                           = (1:N_clusters)';
clusterID                           = repelem(clusterID, T);
year                                = (1:T)';
year                                = 1980+repmat(year,N_clusters,1);

%% define event time now
u1                                  = rand(N_clusters,1);
treated                             = u1>0.2; %80% of cross-sectional units received the policy. 
event_year                          = randi([min(year), max(year)], N_clusters, 1);
treated                             = treated(clusterID); %in the state by year panel;
event_year                          = event_year(clusterID); %in the state by year panel;
event_year                          = treated.*event_year; %0 indicates state was not treated;

%% event-study coefficients (binned at -8 and 8 in this case)
time_rel_event                     = (year-event_year);
time_rel_event(time_rel_event<-wins_value)= -wins_value;
time_rel_event(time_rel_event>wins_value)= wins_value;
lb                                 = -wins_value;
ub                                 = wins_value;
position                           = 1;
D                                  = sparse(size(time_rel_event,1),wins_value*2+1);
for k=min(time_rel_event):max(time_rel_event)
        sel                        = time_rel_event == k;
	    D(sel,position)            = 1;
	    position                   = position + 1;
end

%% create additional controls
rng(3815) %%fix seed to approximate the fixed X design considered in KSS
X                                  = rand(size(treated,1),N_X);
coeff_X                            = [rand(3,1); zeros(N_X-3,1)]; %only three coefficients diff from zero

[~,~,X_T]                          = unique(year,'rows','stable');
X_T 	                           = sparse((1:size(treated,1))',X_T',1,size(treated,1),max(X_T));
coeff_time                         = rand(1,1) + (1:T)'*0.03;

[~,~,X_C]                          = unique(clusterID,'rows','stable');
X_C 	                           = sparse((1:size(treated,1))',X_C',1,size(treated,1),max(clusterID));
coeff_cluster                      = randn(N_clusters,1);


%% created event study coeffs
coeff                              = [0.1.*ones(wins_value,1); 0.2.*ones(wins_value+1,1)];%%parallel trends, then a jump (equal to 0.1)
trend                              = [0.*ones(wins_value,1); 0.01*((1:wins_value+1)'-1)];
coeff                              = coeff + trend;
%disp('First Column: Event Time; Second Column: Event Study Coefficients')
%[(1:size(coeff,1))'-1-wins_value coeff-coeff(wins_value)]

%% serially correlated error term
rng(seed_s) %%this is the part that gets redrawn in each MC draw.
e0                                 = randn(N_clusters,1);
AR_coefficient                     = 0.9;
sigma_WN                           = 0.0001;
e                                  = [e0 zeros(N_clusters,T-1)];  % Initial value as a random number
for t = 2:T
    e(:,t)                         = AR_coefficient * e(:,t-1) + sqrt(sigma_WN)*randn(N_clusters,1);
end 
e_long                             = reshape(e.', [], 1);

%% create the outcome
y                                  = D*coeff + X*coeff_X + X_T*coeff_time  + X_C*coeff_cluster + e_long;

%% normalize event-time dummies
D(:,wins_value)                    = [];
%% save output
%out                                = [y,treated,event_year,clusterID,year,X];
%s=['dgp_panel.csv'];
%dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);





end