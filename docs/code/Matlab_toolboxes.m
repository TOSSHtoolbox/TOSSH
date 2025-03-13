%% checks which Matlab toolboxes are used in TOSSH

% you need to be in the directory above TOSSH
files = dir('./TOSSH/TOSSH_code/*/*.m');

i = 0;
for file = files'
    i = i+1;
    nameList{i} = file.name;
    [fList,pList] = matlab.codetools.requiredFilesAndProducts(nameList{i});
    fList_all{i} = fList'; pList_all{i} = pList'; 
end