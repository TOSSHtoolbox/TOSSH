%% Creates code documentation using M2HTML toolbox

% add directories for functions to path 
% (you need to be in the directory above TOSSH)
if exist('./m2html') == 7
    addpath(genpath('./m2html')); % Matlab to HTML 
else
    error('M2HTML toolbox needed. Can be downloaded from https://www.artefact.tk/software/matlab/m2html/ and should be in a folder named m2html in the same directory.')
end

% check which Matlab version and which toolboxes are installed
ver

% create HTML documentation
m2html('mfiles','TOSSH/TOSSH_code/',...
    'htmldir','TOSSH/docs/matlab/TOSSH_code',...
    'recursive','on')

m2html('mfiles','TOSSH/example/',...
    'htmldir','TOSSH/docs/matlab/example',...
    'recursive','on')
