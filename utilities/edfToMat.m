


function [] = edfToMat(directory,saveFolder)

if ieNotDefined('saveFolder'); saveFolder = directory; end
s = dir(directory);
s2 = s(arrayfun(@(x) ~strcmp(x.name(1),'.'),s));
arrayStringS2  = {s2(:).name};
stimfiles = arrayStringS2(endsWith(arrayStringS2,'.mat'));

for iRun = 1:length(stimfiles)
    cd(directory)
    
    e = getTaskEyeTraces5(stimfiles{iRun});
    edf = e.edf;
    cd(saveFolder)
    
    %idx = find(ismember(e.eyeTrackerFilename,'19')); % CSB, uncomment if it breaks your version
    fileName = e.eyeTrackerFilename;%(idx(1):idx(1)+7); % CSB, uncomment if it breaks your version
    fileName = strcat(fileName, '.mat');
    save(fileName, 'e');
end
end