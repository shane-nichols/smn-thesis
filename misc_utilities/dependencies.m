function depFiles = dependencies(files)
% returns a cell-array of all dependent files
% input files is a single filename or a cell-array of filenames
%  depFiles = dependencies({'myFun1.m' , 'myFun2.m'})
depFiles = matlab.codetools.requiredFilesAndProducts(files);
depFiles = depFiles.';