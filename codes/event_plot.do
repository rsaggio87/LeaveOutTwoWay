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
		
			keep if coeff_ut !=.
			keep coeff_ut lower_ut upper_ut se_ut x 
			drop if x == -1
			gen rows 		= _n
			save "$fileRESULTS",replace
			
*			twoway (scatter coeff_ut x, mcolor(dknavy) lcolor(dknavy) msymbol(square) connect(direct)) ///
*			(rcap lower_ut upper_ut x, lcolor(dknavy)), ///
*			xlabel($lb(1)$ub) legend(off) xtitle("Year Relative to Event") ytitle("Effect and 95% Confidence Intervals") xline(-0.5, lcolor(red) lpattern(dash))  graphregion(color(white)) bgcolor(white)
