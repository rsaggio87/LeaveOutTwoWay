function MakeHierarchyMex

% create mex files
mextocreate = {'adjacency_cmg','diagconjugate','forest_components','graphprofile','laplacian2','perturbtril','splitforest','update_groups','vpack'};
for k=1:length(mextocreate)
    mexstr = sprintf('mex -g -largeArrayDims ../../../Source/Hierarchy/%s.c ',mextocreate{k});
    eval(mexstr);
end

