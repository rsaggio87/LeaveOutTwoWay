
addpath(pwd); 
cd MATLAB/Hierarchy/MexFunctions;
MakeHierarchyMex;
addpath(pwd);
cd ..;
addpath(pwd); 
cd ../..; 
cd MATLAB/Solver;
MakeSolverMex;
addpath(pwd);
cd ../..;
cd MATLAB;
addpath(pwd); 
cd ..;
savepath;

