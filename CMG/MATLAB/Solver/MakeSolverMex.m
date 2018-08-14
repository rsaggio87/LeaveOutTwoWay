function MakeSolverMex

if ispc
    sep = '\';
    obj = 'obj';
else
    sep = '/';
    obj = 'o';
end

include = sprintf('-I..%s..%sInclude',sep,sep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compiling in single precision 
% create object files
objstr=' ';
objtocreate = {'vpv','vmv','vpvmv','vvmul','rmvec','ldl_solve','sspmv','preconditioner'};
for k=1:length(objtocreate)
    filestr = sprintf('..%s..%sSource%sSolver%s%s.c',sep,sep,sep,sep,objtocreate{k});
    mexstr = sprintf('mex -largeArrayDims -DSINGLE_PR %s -c %s',include,filestr);
    objstr = sprintf('%s %s.%s',objstr,objtocreate{k},obj);
    eval(mexstr);
end

% create mex files
mextocreate = {'mx_s_preconditioner'};
for k=1:length(mextocreate)
    mexstr = sprintf('mex -largeArrayDims -DSINGLE_PR %s %s.c %s',include,mextocreate{k},objstr);
    eval(mexstr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compiling in double precision 
% create object files
objstr=' ';
objtocreate = {'vpv','vmv','vpvmv','vvmul','rmvec','ldl_solve','sspmv','preconditioner'};
for k=1:length(objtocreate)
    filestr = sprintf('..%s..%sSource%sSolver%s%s.c',sep,sep,sep,sep,objtocreate{k});
    mexstr = sprintf('mex -largeArrayDims %s -c %s',include,filestr);
    objstr = sprintf('%s %s.%s',objstr,objtocreate{k},obj);
    eval(mexstr);
end

% create mex files
mextocreate = {'mx_d_preconditioner'};
for k=1:length(mextocreate)
    mexstr = sprintf('mex -largeArrayDims %s %s.c %s',include,mextocreate{k},objstr);
    eval(mexstr);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up
delstr = sprintf('delete  %s',objstr);
eval(delstr);