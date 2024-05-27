cd "/Users/raffaelesaggio/Documents/GitHub/LeaveOutTwoWay/data"
use "py_final_1985_2001_veneto_only.dta",clear

keep if year>=1996 & year<=2001
keep if prov == "RO"

keep id firmid year log_dailywages  age prov
gen u1 = runiform()
gen union_status = 0
replace union_status = 1 if u1 >=0.70

rename id id_orig

egen id = group(id_orig union_status)



export delimited id firmid id_orig year log_dailywages union_status age  using "unions_SAKM", replace novarnames nolabel
