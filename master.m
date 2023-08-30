%% paths
clc
clear
restoredefaultpath
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay' %% update this to path where you saved the Github package
path(path,'codes'); %this contains the main LeaveOut Routines.

%% parallel envinr (you can uncommment all of this)
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %TABLE 2: Variance Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
%% read data
namesrc='data/akm_pe_fe_kss_ready.csv'; %path to data
data=importdata(namesrc); %import data
id=data(:,1); %worker identifiers
firmid=data(:,2); %firm identifiers (unrestricted interaction of firmid and union dummy)
firmid_orig=data(:,3); %original firm identifiers
year=data(:,4);
y=data(:,5); % outcome variable
union_status = data(:,6); %0 for non-union job at firm where union share is 0%; 1 for union job; 2 for non-union job at firm with non-zero share of union jobs
age=data(:,7);
clear data

%% age polynomials
age2 = ((age-40)/40).^2;
age3 = ((age-40)/40).^3;

%% create year effects 
[~,~,controls] = unique(year,'rows','stable');
controls 	   = sparse((1:size(y,1))',controls',1,size(y,1),max(controls));			
controls(:,end-1)=[]; %omit one year
controls       = [controls age2 age3];

%% options for KSS
type_of_algorithm='JLA'; 
leave_out_level='matches';
simulations_JLA=100;
call_results='results__KSS_unions';

%% compute table
for type_decomposition=0:2 
tabella(:,type_decomposition+1) = leave_out_KSS(y,id,firmid,firmid_orig,controls,leave_out_level,type_of_algorithm,simulations_JLA,call_results,union_status,type_decomposition); %%calling this will print the variance decomposition for non-union jobs.
end
tabella=[tabella(:,2) tabella(:,3) tabella(:,1)];
tabella=full(tabella);
s=['tables/VARIANCE_DECOMPOSITION.csv'];
dlmwrite(s, tabella, 'delimiter', '\t', 'precision', 16);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %TABLE 3: Oaxaca Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read src files
csv_with_Z                                                = 'data/net_surplus_no_public.csv'; %this must be .csv file where a row correspond to a firm, column 1 reports firmid_orig, column 2 its net surplus, column 3 its industry dummy.here as baseline i am considering a file that excludes firms in the public sector.
csv_with_leave_out_connected_set                          = 'data/results__KSS_unions.csv'; %don't change this.
label_results                                             = 'oaxaca_TABLE_no_public'; %use this so that we can discriminate results that include or not include the public sector. just change path in line 54 and line 56 to have results that include the public sector.

for control_for_industry=0:1
    for type_oaxaca_decomposition=1:3                          
    tabella                                               = oaxaca_projection(csv_with_Z,csv_with_leave_out_connected_set,control_for_industry,type_oaxaca_decomposition);
    net_surplus                                           = tabella(5,:)';
    rent_sharing                                          = tabella(10,:)';
    tabella(13,1)                                         = rent_sharing(1)*net_surplus(2); %Bargaining component (Rent-Sharing coeff of non-union x Surplus of Union)
    tabella(13,2)                                         = rent_sharing(2)*net_surplus(2); %Bargaining component (Rent-Sharing coeff of union x Surplus of Union)
    tabella(14,1)                                         = rent_sharing(1)*net_surplus(1); %Sorting component (Rent-Sharing coeff of non-union x Surplus of Non-Union)
    tabella(14,2)                                         = rent_sharing(1)*net_surplus(2); %Sorting component (Rent-Sharing coeff of non-union x Surplus of Union)
    tabella                                               = [tabella tabella(:,2)-tabella(:,1)];
    SE_cons                                               = tabella(9,:)';
    tabella(9,3)                                          = sqrt(SE_cons(1)^2+SE_cons(2)^2); %%fixing the SEs on the table, no covariance under assumption of error term being correlated only within jobs.
    SE_rs                                                 = tabella(11,:)';
    tabella(11,3)                                         = sqrt(SE_rs(1)^2+SE_rs(2)^2);

    if control_for_industry == 0
        nome_risultati                                    = ['tables/' label_results 'no_controls_for_industry' 'type_decomposition_____' num2str(type_oaxaca_decomposition) '.csv'];
    end

    if control_for_industry == 1
        nome_risultati                                    = ['tables/' label_results 'controls_for_industry'    'type_decomposition_____' num2str(type_oaxaca_decomposition) '.csv'];
    end

    %% print table
    dlmwrite(nome_risultati, tabella, 'delimiter', '\t', 'precision', 16);    
    end
end


