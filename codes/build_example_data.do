********************************************************************************

*		LOAD UP DATA

********************************************************************************
cd "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay/"
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear

l in 1/5

/*

Application is the effects of Medicaid expansions on insurance coverage using publicly-available data from the ACS. This analysis is similar in spirit to that in Carey, Miller, and Wherry (2020), although they use confidential data, see also honest DiD Github package by Rambachan and Roth (2022) which also uses this as their leading application.

The data is a state-level panel with information on health insurance coverage and Medicaid expansion. The variable dins shows the share of low-income childless adults with health insurance in the state. The variable yexp2 gives the year that a state expanded Medicaid coverage under the Affordable Care Act, and is missing if the state never expanded.
*/
********************************************************************************

*		SET UP THE EVENT STUDIES

********************************************************************************
sum year yexp2
gen treated = 1 
replace treated = 0 if yexp2==.

*inspect year of event
rename yexp2 event_year
tab event_year
gen time_rel_event 		= year-event_year
sum time_rel_event

global lb = -6 // winsorize pre-event coefficients at -6
global ub = 4  // winsorize post-event coefficients at 4


replace time_rel_event	=$ub if time_rel_event>=$ub   & time_rel_event!=.
replace time_rel_event	=$lb if time_rel_event<=$lb   & time_rel_event!=.

forval k=$lb/$ub {
		local auxname=`k'-$lb
		
		gen 	D`auxname'	  = 0
		replace D`auxname'	  = 1 	  if time_rel_event==`k' & treated == 1
		
		distinct stfips  time_rel_event if time_rel_event==`k' // tell me how many states for a specific value of time_rel_event
}

local norma = -$lb - 1
********************************************************************************

*		EXPORT DATASET TO MATLAB		

********************************************************************************
keep dins D* stfips year
order dins D* stfips year 
drop D`norma' // otherwise MATLAB code would crash because of collinearity issue. 
export delimited using "data_MEDICAID.csv", replace novarnames nolabel
reghdfe dins D*, absorb(stfips year) cluster(stfips) noconstant
global fileRESULTS = "results_MEDICAID"
do "codes/event_plot" // compute event stdu plot and save coefficients


