function tabella=oaxaca_projection(csv_with_Z,csv_with_leave_out_connected_set,control_for_industry,type_oaxaca_decomposition)

%% load the file with information on net surplus and industry dummies
namesrc=[csv_with_Z]; %path to data
data=importdata(namesrc); %import data
firmid_orig=data(:,1); %original firm identifiers
net_surplus=data(:,2); %net surplus
ind_dummies=data(:,3); %industry dummies
elist1= [firmid_orig net_surplus ind_dummies];%% Original Firmid; N of Union Jobs; Id for AKM
clear data

%% load the personfile saved by "leave_out_KSS"
namesrc=[csv_with_leave_out_connected_set]; %path to data
data=importdata(namesrc); %import data
N = max(data(:,8));
J = max(data(:,9));

%% merge information
LIST_BASE 				   = array2table(elist1,...
    								'VariableNames',{'firmid_original','net_surplus','ind_dummies'});
LIST_SEL 				   = array2table(data,...
        							'VariableNames',{'id_old','firmid_old','firmid_original','union_status','y','fe','pe','id','firmid'});
merge					   = outerjoin(LIST_BASE,LIST_SEL);
merge					   = outerjoin(LIST_BASE,LIST_SEL);
merge 					   = table2array(merge);
data					   = merge(~any(isnan(merge),2),:); %% keep if _merge ==3 --->retain only person-year obs associated with firms present in 'csv_with_Z'
clear merge

%% extract the Z variable for lincom + additional variables that we will eventually summarize
net_surplus                = data(:,2);
ind_dummies                = data(:,3);
union_status               = data(:,3+4);
y                          = data(:,3+5);
fe                         = data(:,3+6);
pe                         = data(:,3+7);
firmid                     = data(:,3+9);

%% set up things to run lincom
load('data/designMatrix results__KSS_unions')
sigma_i                   = myMatrix(:,1);
yall                      = myMatrix(:,2); %includes the outcome for all observations, not just those with, say, union_status==1, or union_status==0;
X                         = myMatrix(:,3:end);

%% now tell me which type of decomposition you are doing

if type_oaxaca_decomposition==1
    treatment              = union_status;
end

if type_oaxaca_decomposition==2
    treatment              = union_status;
    treatment(treatment==2)= 0;
end

if type_oaxaca_decomposition==3
    treatment              = union_status;
    treatment(treatment==2)= 1;
end

%% create industry dummies 
if control_for_industry == 1
    [~,~,ind_dummies]       = unique(ind_dummies);
    ind_dummies 	        = sparse((1:size(ind_dummies,1))',ind_dummies',1,size(ind_dummies,1),max(ind_dummies));	
    ind_dummies(:,end-1)    =[]; %omit one  
end
if control_for_industry == 0
    ind_dummies             =[]; 
end

%% now run lincom separately for different values of the union status variable
tabella                     = zeros(14,2);

for val_treatment=0:1
    sel                     = treatment == val_treatment;
    Nsel                    = sum(sel);
    Transform               = sparse((1:Nsel)',firmid(sel),1,Nsel,J);
    S                       = speye(J-1);
    S                       = [S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
    Transform               = Transform*S;
    Transform               = [sparse(Nsel,N) Transform];
    
    if control_for_industry == 1
        Z                   = [net_surplus(sel) ind_dummies(sel,:)];
    end
    
    if control_for_industry == 0
        Z                   = [net_surplus(sel)];
    end
    
    [~, linear_combination, SE_linear_combination_KSS] = lincom_KSS(yall,X,Z,Transform,sigma_i);

    tabella(1,val_treatment+1)=mean(y(sel));
    tabella(2,val_treatment+1)=-9;
    tabella(3,val_treatment+1)=mean(pe(sel));
    tabella(4,val_treatment+1)=mean(fe(sel));
    tabella(5,val_treatment+1)=mean(net_surplus(sel));
    tabella(6,val_treatment+1)=-9;
    tabella(7,val_treatment+1)=-9;
    tabella(8,val_treatment+1)=linear_combination(1);
    tabella(9,val_treatment+1)=SE_linear_combination_KSS(1);
    tabella(10,val_treatment+1)=linear_combination(2);
    tabella(11,val_treatment+1)=SE_linear_combination_KSS(2);
    tabella(12,val_treatment+1)=-9;
    tabella(13,val_treatment+1)=-9;
    tabella(14,val_treatment+1)=-9;
end
end