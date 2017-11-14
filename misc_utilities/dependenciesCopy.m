function dependenciesCopy(filenames,fullPath)
% same as function 'dependencies' except that here all dependent files are
% copied to the location of fullPath
files = dependencies(filenames);
if exist(fullPath,'dir') == 0
    mkdir(fullPath);
end
for i=1:length(files)
    [~,name,ext] = fileparts(files{i});
    newFile = fullfile(fullPath,[name,ext]);
    copyfile(files{i},newFile)
end
end