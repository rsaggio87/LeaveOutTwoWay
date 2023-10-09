# Leave-Out Cluster Standard Errors

This readme shows how to compute an estimate of the samplig variability of linear regression coefficients that is unbiased even in the presence of serial-correlation of the error term within cluster. This might be particularly useful when estimating a regression model with very few clusters. This [document](https://www.dropbox.com/scl/fi/vxyss0tf3h50lpwrp80c0/metrics.pdf?rlkey=ne9yiquzcj3k9d4itx4vzzlm1&dl=0) describes and provides intuition for the econometric formula used in this package and derived in Remark 3 of [Kline-Saggio-Sølvsten (2020)](https://eml.berkeley.edu/~pkline/papers/KSS2020.pdf)--KSS henceforth.

# The `KSS_SE` Function

The function `KSS_SE` is a `MATLAB` function that takes as input the following:
* `y` outcome variable. Dimension is $N\times 1$.
* `D` treatment(s) of interest. Dimension is $N\times N_{D}$.
* `clusterID` variable that indexes clusters. Dimension is $N\times 1$.
* `indexID` variable that indexes observations within clusters. Dimension is $N\times 1$.
* `controls` additional controls. Dimension is $N\times N_{P}$.

The function computes the KSS leave-out standard errors on the regression coefficients associated with `D` after controlling for `controls`, `clusterID` fixed effects as well as `indexID` fixed effects. These SEs are clustered at the level indexed by `clusterID`.  The inputs  `indexID` and `controls` are optional and can be supplied as empty arrays `[]`.

We now demonstrate the functioning of  `KSS_SE` in the context where one is interested in fitting an event study model of the form

$$y_{it} = \alpha_{i} + \lambda_{t} + \sum_{k=a}^{b}D_{it}^{k}\theta_{k}+X_{it}'\gamma + r_{it}$$
where $\alpha_{i}$ are, say, state fixed effects; $\lambda_{t}$ are year fixed effects; $D_{it}^{k}$ are event study indicators and $X_{it}$ are some time-varying controls. 

# Building and Exporting the Data To Matlab

We first load up ACS public available data on health insurance using `Stata`. 


```matlab
import stata_setup
stata_setup.config("/Applications/STATA", "se")
```

    
      ___  ____  ____  ____  ____ ®
     /__    /   ____/   /   ____/      18.0
    ___/   /   /___/   /   /___/       SE—Standard Edition
    
     Statistics and Data Science       Copyright 1985-2023 StataCorp LLC
                                       StataCorp
                                       4905 Lakeway Drive
                                       College Station, Texas 77845 USA
                                       800-STATA-PC        https://www.stata.com
                                       979-696-4600        stata@stata.com
    
    Stata license: Unlimited-user network, expiring 19 Aug 2024
    Serial number: 401809301518
      Licensed to: Raffaele Saggio
                   UBC
    
    Notes:
          1. Unicode is supported; see help unicode_advice.
          2. Maximum number of variables is set to 5,000 but can be increased;
              see help set_maxvar.



```matlab
%%stata
cd "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay"
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
l in 1/5
tab year
tab yexp2
gen treated = 1 
replace treated = 0 if yexp2==.
```

    
    . cd "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay"
    /Users/raffaelesaggio/Dropbox/LeaveOutTwoWay
    
    . local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
    
    . use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
    
    . l in 1/5
    
         +--------------------------------------------+
         |  stfips   year       dins   yexp2        W |
         |--------------------------------------------|
      1. | alabama   2008   .6814122       .   613156 |
      2. | alabama   2009   .6580621       .   613156 |
      3. | alabama   2010   .6313651       .   613156 |
      4. | alabama   2011   .6563886       .   613156 |
      5. | alabama   2012   .6708115       .   613156 |
         +--------------------------------------------+
    
    . tab year
    
     Census/ACS |
    survey year |      Freq.     Percent        Cum.
    ------------+-----------------------------------
           2008 |         46        8.33        8.33
           2009 |         46        8.33       16.67
           2010 |         46        8.33       25.00
           2011 |         46        8.33       33.33
           2012 |         46        8.33       41.67
           2013 |         46        8.33       50.00
           2014 |         46        8.33       58.33
           2015 |         46        8.33       66.67
           2016 |         46        8.33       75.00
           2017 |         46        8.33       83.33
           2018 |         46        8.33       91.67
           2019 |         46        8.33      100.00
    ------------+-----------------------------------
          Total |        552      100.00
    
    . tab yexp2
    
        Year of |
       Medicaid |
      Expansion |      Freq.     Percent        Cum.
    ------------+-----------------------------------
           2014 |        264       73.33       73.33
           2015 |         36       10.00       83.33
           2016 |         24        6.67       90.00
           2017 |         12        3.33       93.33
           2019 |         24        6.67      100.00
    ------------+-----------------------------------
          Total |        360      100.00
    
    . gen treated = 1 
    
    . replace treated = 0 if yexp2==.
    (192 real changes made)
    
    . 


This is a state-year panel where the variable `dins` is the outcome of interest and `yexp2` measures the year in which Medicaid was expanded in a given state (it is missing for states that did not expand, like Alabama). Note that the panel runs from 2008 and 2019  and most states expanded in 2014.

We now export to matlab `dins` `stfips` `year` and the set of event study indicators in a .csv called `data_MEDICAID.csv`. To do that, we winsorize the event-study indicators at -6 and +4, export the resulting data and save the results after fitting the event-study specification


```matlab
%%stata
rename yexp2 event_year
gen time_rel_event 		= year-event_year
sum time_rel_event
global lb = -6 // winsorize pre-event coefficients at -6
global ub = 4  // winsorize post-event coefficients at 4
replace time_rel_event	=$ub if time_rel_event>=$ub   & time_rel_event!=.
replace time_rel_event	=$lb if time_rel_event<=$lb   & time_rel_event!=.
qui forval k=$lb/$ub {
		local auxname=`k'-$lb
		gen 	D`auxname'	  = 0
		replace D`auxname'	  = 1 if time_rel_event==`k' & treated == 1
}
local norma = -$lb - 1
replace D`norma'=0 // normalize relative to year before implemention
keep dins D* stfips year
order dins D* stfips year 
drop D`norma' // do not import one event study indicator otherwise MATLAB code would crash because of collinearity issue. 
export delimited using "data/data_MEDICAID.csv", replace novarnames nolabel
reghdfe dins D*, absorb(stfips year) cluster(stfips) noconstant
```

    
    . rename yexp2 event_year
    
    . gen time_rel_event              = year-event_year
    (192 missing values generated)
    
    . sum time_rel_event
    
        Variable |        Obs        Mean    Std. dev.       Min        Max
    -------------+---------------------------------------------------------
    time_rel_e~t |        360   -1.166667    3.720754        -11          5
    
    . global lb = -6 // winsorize pre-event coefficients at -6
    
    . global ub = 4  // winsorize post-event coefficients at 4
    
    . replace time_rel_event  =$ub if time_rel_event>=$ub   & time_rel_event!=.
    (22 real changes made)
    
    . replace time_rel_event  =$lb if time_rel_event<=$lb   & time_rel_event!=.
    (20 real changes made)
    
    . qui forval k=$lb/$ub {
    
    . local norma = -$lb - 1
    
    . replace D`norma'=0 // normalize relative to year before implemention
    (30 real changes made)
    
    . keep dins D* stfips year
    
    . order dins D* stfips year 
    
    . drop D`norma' // do not import one event study indicator otherwise MATLAB cod
    > e would crash because of collinearity issue. 
    
    . export delimited using "data/data_MEDICAID.csv", replace novarnames nolabel
    file data/data_MEDICAID.csv saved
    
    . reghdfe dins D*, absorb(stfips year) cluster(stfips) noconstant
    (MWFE estimator converged in 2 iterations)
    
    HDFE Linear regression                            Number of obs   =        552
    Absorbing 2 HDFE groups                           F(  10,     45) =      16.52
    Statistics robust to heteroskedasticity           Prob > F        =     0.0000
                                                      R-squared       =     0.9482
                                                      Adj R-squared   =     0.9411
                                                      Within R-sq.    =     0.4631
    Number of clusters (stfips)  =         46         Root MSE        =     0.0218
    
                                    (Std. err. adjusted for 46 clusters in stfips)
    ------------------------------------------------------------------------------
                 |               Robust
            dins | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
    -------------+----------------------------------------------------------------
              D0 |  -.0034914   .0082004    -0.43   0.672    -.0200079    .0130251
              D1 |  -.0108234   .0069703    -1.55   0.127    -.0248623    .0032156
              D2 |   -.008634   .0057581    -1.50   0.141    -.0202314    .0029633
              D3 |  -.0049002   .0053501    -0.92   0.365    -.0156759    .0058755
              D4 |  -.0065987   .0038092    -1.73   0.090    -.0142707    .0010734
              D6 |   .0445785    .005667     7.87   0.000     .0331645    .0559925
              D7 |   .0643379   .0072982     8.82   0.000     .0496385    .0790372
              D8 |   .0779544   .0082594     9.44   0.000     .0613191    .0945896
              D9 |   .0737635   .0095947     7.69   0.000     .0544388    .0930882
             D10 |   .0752858   .0108269     6.95   0.000     .0534794    .0970922
    ------------------------------------------------------------------------------
    
    Absorbed degrees of freedom:
    -----------------------------------------------------+
     Absorbed FE | Categories  - Redundant  = Num. Coefs |
    -------------+---------------------------------------|
          stfips |        46          46           0    *|
            year |        12           0          12     |
    -----------------------------------------------------+
    * = FE nested within cluster; treated as redundant for DoF computation
    
    . 


## Running KSS_SE in MATLAB
We now switch to MATLAB and begin by setting the path and the parfor environment.


```matlab
%% setting path and parfor env
clc
clear
cd '/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay' %% update this to path where you saved the Github package
path(path,'codes'); %this contains the main LeaveOut Routines.

%% parallel envinr (you can uncommment all of this)
%delete(gcp("nocreate")) %clear parallel envir.
%c = parcluster('local'); %tell me # of available cores
%nw = c.NumWorkers; %tell me # of available cores
%pool=parpool(nw,'IdleTimeout', Inf); %all cores will be assigned to Matlab
```

We now import the data exported by STATA


```matlab
namesrc='data/data_MEDICAID.csv'; %see build_example_data
data=importdata(namesrc); %import data
y=data(:,1); %outcome variable
D=data(:,2:11); %treatment of interest. in the example D corresponds to 9 event_study dummies. note that event dummy for the year before medicaid was implemented is not imported to avoid multicollinearity issues. 
clusterID=data(:,12); %this is a variable that contains the identifier of the cluster. this is the level at which we are going to cluster the standard errors
year=data(:,13); %other dimension of the panel, in this case this is year (data is state by year panel).
controls=[]; %in this example no controls but the code allows the inclusion of controls
clear data

%%next step is just to assign a label to each column of D, this is not
%%necessary.
for k = -6:-2
    labels{k + 7}               = ['Event-Study Coefficient at ' num2str(k)]; % Adjust the index
end

for k = 0:4
    labels{k + 6}               = ['Event-Study Coefficient at ' num2str(k)]; % Adjust the index
end
```

We are ready to launch `KSS_SE` by simply calling


```matlab
%% Run KSS!
[coeff, SE]   = KSS_SE(y,D,clusterID,year,controls,labels); %% this is equivalent in STATA to run reghdfe y D* controls, abs(clusterID year) cluster(clusterID)

%coeff are the event-study coefficients in this application.
%SE are the KSS standard errors for coeff computed by leaving-out cluster c out, as described in Remark 3 of KSS.

out                             = [coeff,SE]; 
s                               = ['data/results_MEDICAID_MATLAB.csv'];
dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); %% saving results in a  .csv
```

    Coefficient on Event-Study Coefficient at -6:  -0.0034914
    Leave-Out Standard Error: 0.0089091
    ******************************************
    Coefficient on Event-Study Coefficient at -5:  -0.010823
    Leave-Out Standard Error: 0.0068348
    ******************************************
    Coefficient on Event-Study Coefficient at -4:  -0.008634
    Leave-Out Standard Error: 0.0075653
    ******************************************
    Coefficient on Event-Study Coefficient at -3:  -0.0049002
    Leave-Out Standard Error: 0.0049889
    ******************************************
    Coefficient on Event-Study Coefficient at -2:  -0.0065987
    Leave-Out Standard Error: 0.0036725
    ******************************************
    Coefficient on Event-Study Coefficient at 0:  0.044579
    Leave-Out Standard Error: 0.0061763
    ******************************************
    Coefficient on Event-Study Coefficient at 1:  0.064338
    Leave-Out Standard Error: 0.0065041
    ******************************************
    Coefficient on Event-Study Coefficient at 2:  0.077954
    Leave-Out Standard Error: 0.008109
    ******************************************
    Coefficient on Event-Study Coefficient at 3:  0.073763
    Leave-Out Standard Error: 0.010457
    ******************************************
    Coefficient on Event-Study Coefficient at 4:  0.075286
    Leave-Out Standard Error: 0.011967
    ******************************************


# Interpreting the Output

The code `KSS_SE` prints the event-study coefficients (normalized relative to the year expansion of Medicaid) along with KSS standard errors, clustered at the state level. As shown in KSS, the (square) of the standard errors printed by `KSS_SE`, represent an unbiased estimate of true sampling variable of these OLS estimates. The so-called "White" standard errors that are reported by standard packages (e.g. `reghdfe`)  are consistent but are typically biased in finite samples, especially when the number of clusters is small. We demonstrate this point by means of a montecarlo simulation below.

# Montecarlo Exercise

Consider a simple DGP of the type 

$$y_{it} = \alpha_{i} + \lambda_{t} + \sum_{k=a}^{b}D_{it}^{k}\theta_{k}+X_{it}'\gamma + r_{it}$$

where $r_{it}$ is drawn from an AR(1) process with AR(1) coefficient equal to 0.9 and normal white noise errors. In this DGP, there are 30 clusters , 20 time periods, 10 extra controls and a total of 30 regression coefficients of interest (i.e $a$=-15 and $b$=15 in the event-study regression above). We simulate from this DGP and print the sampling variability of $\theta_{k}$ across these 1000 simulations along with the estimate of this sampling variability based on the leave-out formula of KSS as well as estimates based on the traditional "White" cluster-robust standard errors that are typically printed by standard packages such as `reghdfe`. 


```matlab
%% some parameters for simulated data
wins_value                      = 15; %event studies binned at -wins_value and wins_value
N_clusters                      = 30; %number of clusters.
N_X                             = 10; %number of extra controls.
T                               = 20; %number of time periods.
S                               = 1000; %number of montercarlo draws
coeff_sim                       = zeros(S,wins_value*2); %this will contain the estimated \theta_k for a given MC draw 
V_sim                           = zeros(S,wins_value*2); %this will contain the estimated sampling variability of \theta_k based on KSS
V_sim_naive                     = zeros(S,wins_value*2); %this will contain the estimated sampling variability of \theta_k based on White 
parfor s = 1:S
[y, D, clusterID, year, X ]     = create_DGP(wins_value,N_clusters,N_X,T,s); %this takes a draw from the DGP described above
[coeff_sim(s,:),SE,SE_naive]    = KSS_SE(y,D,clusterID,year,X); %run KSS_SE on the simulated data.
V_sim(s,:)                      = SE.^2; %%store the squared of the SE based on KSS
V_sim_naive(s,:)                = SE_naive.^2; %%store the squared of the SE based on White
end
true_coeff                      = [0.1.*ones(wins_value,1); 0.2.*ones(wins_value+1,1)];%%parallel trends, then a jump (equal to 0.1)
trend                           = [0.*ones(wins_value,1); 0.01*((1:wins_value+1)'-1)];
true_coeff                      = true_coeff + trend; %% this is how the true event-study coefficients are defined in the DGP, see line 45 of create_DGP.
true_coeff                      = true_coeff-true_coeff(wins_value); %%OLS is normalized relative to year before event.
true_coeff(wins_value)          = []; %omit year before event
C                               = cov(coeff_sim); %%VCM of OLS estimates across MC draws
true_V                          = diag(C); %%get variances
coeff_sim                       = mean(coeff_sim)'; %mean of OLS
V_sim                           = mean(V_sim)'; %mean of leave-out estimates of the (squared) SE
V_sim_naive                     = mean(V_sim_naive)'; %mean of White estimates of the (squared) SE
myArray                         = [true_coeff coeff_sim  true_V V_sim V_sim_naive];
a2t                             = array2table(myArray,"VariableNames",["True Event-Study Coefficients","Mean of OLS","Variance of OLS" "Leave-Out Estimate of Variance of OLS" "White Estimate of Variance of OLS"]);
disp(a2t)
```

        True Event-Study Coefficients    Mean of OLS    Variance of OLS    Leave-Out Estimate of Variance of OLS    White Estimate of Variance of OLS
        _____________________________    ___________    _______________    _____________________________________    _________________________________
                       0                   0.020049          0.37206                       0.35382                                0.26345            
                       0                    0.02895          0.23803                       0.22334                                0.17832            
                       0                    0.02593          0.18583                       0.17773                                0.14673            
                       0                   0.020056          0.14956                       0.14327                                0.12008            
                       0                   0.017036          0.11813                       0.11386                               0.095974            
                       0                    0.01336         0.095366                      0.090648                               0.076853            
                       0                   0.011525         0.075004                      0.071976                               0.060818            
                       0                  0.0097184         0.057437                       0.05546                               0.046774            
                       0                  0.0095178         0.042274                      0.041374                               0.035192            
                       0                  0.0064399         0.030653                      0.029811                               0.025314            
                       0                  0.0027579         0.020411                       0.02049                               0.017475            
                       0                  0.0029581         0.012726                      0.013066                               0.011109            
                       0                 0.00092733        0.0067026                     0.0072042                              0.0060794            
                       0                 -0.0014629        0.0026599                     0.0029761                              0.0025117            
                     0.1                   0.096903         0.002644                     0.0028592                              0.0023614            
                    0.11                    0.10577        0.0061617                      0.005905                              0.0049026            
                    0.12                    0.11696         0.011424                      0.010805                              0.0089748            
                    0.13                    0.12547         0.018587                      0.017333                               0.014444            
                    0.14                    0.13501         0.028104                      0.025716                               0.021412            
                    0.15                    0.14538         0.039463                      0.035931                               0.029969            
                    0.16                    0.15581         0.052477                      0.048072                               0.040015            
                    0.17                     0.1655         0.068675                      0.062081                               0.051882            
                    0.18                     0.1746          0.08745                      0.077587                               0.065276            
                    0.19                    0.18554          0.10709                      0.095506                               0.080399            
                     0.2                    0.19529           0.1307                       0.11711                               0.097787            
                    0.21                    0.20481          0.15522                       0.14028                                0.11662            
                    0.22                    0.21794          0.18407                        0.1665                                0.13744            
                    0.23                    0.22939          0.21667                       0.19511                                0.15997            
                    0.24                    0.23935          0.24664                       0.22852                                0.18625            
                    0.25                    0.24774          0.33507                       0.30833                                 0.2397            


One can see from the output above that the leave-out estimate of the sampling variability are approximately unbiased while the White SEs are downward biased. This is particularly evident when looking at the tails of the event-study coefficients (e.g last two-three rows of the table above). This is to be expected given the unbalanced nature of the design as there are very few clusters for which we are able to observe outcomes, say, 15 years after the implementation of the policy. 
