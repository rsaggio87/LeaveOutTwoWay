********************************************************************************

*		LOAD UP DATA

********************************************************************************
cd "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay/"
import delimited "dgp_panel.csv",clear
local wins_value = 5 // set bins
********************************************************************************

*		GET READY TO RUN EVENT-STUDY

********************************************************************************
rename v1 y
rename v2 treated
rename v3 event_year
rename v4 cluster_id
rename v5 year

gen time_rel_event 		= year-event_year
replace time_rel_event	=`wins_value' if time_rel_event>=`wins_value'   & time_rel_event!=.
replace time_rel_event	=-`wins_value' if time_rel_event<=-`wins_value' & time_rel_event!=.

global lb = -`wins_value'
global ub = `wins_value'


forval k=$lb/$ub {
		local auxname=`k'-$lb
		
		gen 	D`auxname'	  = 0
		replace D`auxname'	  = 1 	  if time_rel_event==`k' & treated == 1

}
local norma = `wins_value'-1
replace D`norma'=0

********************************************************************************

*		RUN EVENT-STUDY

********************************************************************************
reghdfe y D* v*, abs(cluster_id year) cluster(cluster_id)
stop
mat coeff1=[_b[D0]]
			mat se=[_se[D0]]
			mat upper=coeff1+((1.96)*se)
			mat lower=coeff1-((1.96)*se)
			mat x=$lb
			local max_bins=-$lb+$ub
		
			forval k=1/`max_bins'{
				mat coeff1=[coeff1\ _b[D`k']]
				mat se=[se\ _se[D`k']]
				mat upper=coeff1+((1.96)*se)
				mat lower=coeff1-((1.96)*se)
				mat x=[x\ $lb+`k']
			}
		
			mat ut_data=[coeff1,upper,lower,se,x]
			svmat ut_data
			rename ut_data1 coeff_ut
			rename ut_data2 upper_ut
			rename ut_data3 lower_ut
			rename ut_data4 se_ut
			rename ut_data5 x			
********************************************************************************

*		COMPUTE GRAPH THAT COMBINES STANDARD ERRORS

********************************************************************************			
keep if coeff_ut !=.
keep coeff_ut lower_ut upper_ut x 
drop if x == -1
gen rows 		= _n
save "coeff_naive",replace
import delimited "results_AFP.csv",clear
gen rows 		= _n
rename v1 coeff_KSS
rename v2 SE_KSS
gen lower_KSS = coeff_KSS-1.96*SE_KSS
gen upper_KSS = coeff_KSS+1.96*SE_KSS
merge 1:1 rows using "coeff_naive", nogen
sort x
twoway (scatter coeff_ut x, mcolor(dknavy) lcolor(dknavy) msymbol(square) connect(direct) legend(label(1 "Cluster White SE"))) ///
	   (rcap lower_ut upper_ut x, lcolor(dknavy) legend(label(2 ""))) /// 
	   (scatter coeff_KSS x, mcolor(cranberry) lcolor(cranberry) msymbol(square) connect(direct) legend(label(3 "Cluster KSS SE"))) ///
	   (rcap lower_KSS upper_KSS x, lcolor(cranberry) legend(label(4 ""))), /// 
		xlabel($lb(1)$ub) xtitle("Year Relative to Event") legend(order(1 3)) legend(ring(0) position(6) rows(2)) xline(-0.5, lcolor(red) lpattern(dash))  graphregion(color(white)) bgcolor(white) note("95% CI shown in the graph")

		
stop

