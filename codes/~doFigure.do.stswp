cd "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay"
import delimited "data/results_MEDICAID_MATLAB.csv",clear
gen rows 		= _n
rename v1 coeff_KSS
rename v2 SE_KSS
gen lower_KSS = coeff_KSS-1.96*SE_KSS
gen upper_KSS = coeff_KSS+1.96*SE_KSS
merge 1:1 rows using "data/results_MEDICAID", nogen
sort x
twoway (scatter coeff_ut x, mcolor(dknavy) lcolor(dknavy) msymbol(square) connect(direct) legend(label(1 "CI based on Cluster-Robust White SE"))) ///
	   (rcap lower_ut upper_ut x, lcolor(dknavy) legend(label(2 ""))) /// 
	   (scatter coeff_KSS x, mcolor(cranberry) lcolor(cranberry) msymbol(square) connect(direct) legend(label(3 "CI based on Cluster-Robust KSS SE"))) ///
	   (rcap lower_KSS upper_KSS x, lcolor(cranberry) legend(label(4 ""))), /// 
		xlabel(-6(1)4) ytitle("Effect on share of low-income childless adults with health insurance", size(small)) xtitle("Year Relative to Expansion of Medicaid") legend(order(1 3)) legend(ring(0) position(11) rows(2)) xline(-0.5, lcolor(red) lpattern(dash))  graphregion(color(white)) bgcolor(white) note("Graph shows the effect of Medicaid expansions on insurance coverage using a state-year panel from ACS" "Red CI based on KSS-SEs clustered at the state-level")

keep coeff_KSS SE_KSS se_ut x
rename coeff_KSS coeff
rename se_ut SE_WHITE
list
save "data/final_TABLE_comparison_KSS",replace

