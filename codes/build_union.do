****load up the original dataset used to test KSS
import delimited "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay/data/lincom.csv",clear
keep v1 v2 v3 v5

rename v1 id
rename v2 firmid
rename v3 year
rename v5 y

egen match_id = group(id firmid)

*randomly make 20% of jobs unionized
preserve
gen uno = 1
collapse uno, by(match_id)
drop uno
generate u1 = runiform()
gen union = 0
replace union = 1 if u1>0.80
drop u1
save "/Users/raffaelesaggio/Dropbox/Downloads/unions_SIM",replace
restore
merge m:1 match_id using "/Users/raffaelesaggio/Dropbox/Downloads/unions_SIM",nogen

*calculate union share at firm level
bys firmid: egen union_share = mean(union)
sum union_share,d


*calculate the firmid x union status 
egen firmid_by_union = group(union firmid)

*create the variable union_status now
gen 		union_status = 1 if union == 1
replace		union_status = 0 if union_share==0
replace		union_status = 2 if union_share>0 & union == 0
tab union_status


*final touches
generate age = runiformint(20, 60)
order id firmid_by_union firmid year y union_status age

xtset id year
export delimited using "/Users/raffaelesaggio/Dropbox/LeaveOutTwoWay/data/unions.csv", novarnames replace
