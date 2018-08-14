A. Installation:

Call MakeCMG.m


*********************************************************
Note: 
MakeCMG adds to the MATLAB path the following directories:

.
./MATLAB/Hierarchy/Mexfunctions;
./MATLAB/Hierarchy;
./MATLAB/solver; 
./MATLAB;

It then saves the path. In shared installations of MATLAB
the path is not saved for next session. So when you re-run
MATLAB make sure that the directories are on the path. 

*********************************************************

B. Runnning CMG:

------------------------------------------------------------------
CAVEATS:

Currently the solver supports only non-positive off-diagonal elements.
When the solver runs into a dense matrix it will issue a warning. 
The development of a sparsification routine is under way. 
------------------------------------------------------------------

A simple scenario. Assume you want to solve many systems
with a matrix A, which is diagonally dominant with 
negative off-diagonals. 

Run one time:

>> pfun = cmg_sdd(A); 

Then for each different b-side, run:

>> x = pcg(A, b, tol,iter, pfun);  % help pcg for documentation


For help:
>> help cmg_sdd




