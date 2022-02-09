function [d] = dataAnalysis(directory)
% dataAnalysis.m
%
%      Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%      Date:    09/26/2020
%
%      Purpose: Estimates trial-averaged pupil response, saccade-locked
%               pupil response, and saccade rate function.
%
%       Input:  "directory", path containing .mat files - both stimfiles and converted edf file
%
%      Outputs: "d", data struct containing pupil trial averages,
%               deconvolved saccade-locked pupil responses, and
%               saccade rate functions for each run (and separately for
%               correct and incorrect trials).
%
%     Todo:     - add a flag for each session
%               - automatically run edfToMat in directory if it hasn't been
%                 run already yet.


if ieNotDefined('directory'); directory = pwd; end

pupilAvg = [];
sacRate = [];
sacIrf =[];
baseline = [];
lookBeforeOnset = 0;

% make this so that it only finds and loads files that aren't HIDDEN- with
% a "." before the filename. Also it only load mat files.
s = dir(directory);
s2 = s(arrayfun(@(x) ~strcmp(x.name(1),'.'),s));
arrayStringS2  = {s2(:).name};
stimfiles = arrayStringS2(endsWith(arrayStringS2,'.mat'));
stimfileInds = contains(stimfiles,'stim');
stimfiles = stimfiles(stimfileInds);
cd(directory)

pupTrace = cell(1,length(stimfiles));
trInds = cell(1,length(stimfiles));
accuracyByTrial = cell(1,length(stimfiles));

for iRun = 1:length(stimfiles)
    
    % load the stimfile
    s = load(stimfiles{iRun});
    
    % get the seglen and compute trial duration
    try
        trialLength = zeros(1,length(s.task{1,1}.seglenPrecompute.seglen.vals));
        for ii = 1:length(s.task{1,1}.seglenPrecompute.seglen.vals) % for case when the timing is jittered and the ITIs are precomputed
            trialLength(ii) = sum(s.task{1,1}.seglenPrecompute.seglen.vals{ii});
        end
        d.maxTrialLengthSec = max(trialLength);
        d.firstSegLength = s.task{1,1}.seglenPrecompute.seglen.vals{1}(1);
    catch
        d.maxTrialLengthSec = sum(s.task{1,1}.seglen);
        d.firstSegLength = s.task{1,1}.seglen(1);
    end
    
    
    %% load and format eye data
    
    % load edf file
    e = getTaskEyeTraces3(stimfiles{iRun});
    edf = e.edf;
    
  
    tiltStim2 = e.randVars.tiltStim2; %s.stimulus.tiltStim2; % for fixedITI data instead, change to s.stimulus.tiltStim2   .   % for cue validity expt, set to "s.stimulus.tiltStim2(1:120);" % for block alt, set to e.randVars.tiltStim2(1:160)
    tiltStim2(tiltStim2==1)=2;
    tiltStim2(tiltStim2==-1)=1;
    trialCorrect = tiltStim2==e.response;
    accuracyByTrial{iRun} = trialCorrect;

    
    sampleRate = e.edf.samplerate;
    
    % extract time stamps of eye tracker, trial-onset, and button press
    iTkrTime = edf.gaze.time;
    trialStart = edf.mgl.time(edf.mgl.segmentNum==1);
    
    % find the index of the start of the run
    [~, startInd] = min(abs(iTkrTime-trialStart(1)));
    startInd = startInd(1);
    
    % find the last segment length
    try
        lastSegLength = s.task{1,1}.seglenPrecompute.seglen.vals{end}(end);
    catch
        lastSegLength = s.task{1,1}.seglen(end);
    end
    
    % find index of the end of the run
    [~, endInd] = min(abs(iTkrTime-(edf.mgl.time(end)+lastSegLength*1000)));
    endInd = endInd(1);
    
    % store the time indicies
    eyeTrackerTime = startInd:endInd;
    
    % extract the appropriate chunk of eye data
    pupil = edf.gaze.pupil(eyeTrackerTime)';
    
    % extract the position data
    xPos = edf.gaze.x(eyeTrackerTime)';
    yPos = edf.gaze.y(eyeTrackerTime)';
    
    % convert to degrees of visual angle
    w = s.myscreen.screenWidth;
    h = s.myscreen.screenHeight;
    xPix2Deg = s.myscreen.imageWidth/w;
    yPix2Deg = s.myscreen.imageHeight/h;
    xPos = (xPos-(w/2))*xPix2Deg;
    yPos = ((h/2)-yPos)*yPix2Deg;
    
    %% detect (micro)saccades
    
    % now look for microsaccades
    eyepos = cat(2, xPos, yPos);
    
    % compute velocity
    vel = vecvel(eyepos, sampleRate, 3);
    
    % detect microsaccades (using method fo Engbert & Mergenthaler 2006)
    mindur = 7;
    vthresh = 8;
    sacs = microsacc(eyepos, vel, vthresh, mindur);
    
    sacTimes{iRun} = sacs; % save out saccade times to fit inter-saccadic interval distribution later
    
    %%
    edf.inds = edf.inds-startInd+1;
    edf.inds(:,2) = edf.inds(:,2) - 1;
    edf.inds(end,2) = length(pupil);
    trInds{iRun} = edf.inds;
    
    %% estimate saccade-locked pupil response (deconvolve pupil size to saccade onset times)
    
    sacTime = sacs(:,1); % determine the indexes that saccades begin
    downsampleRate = 5;
    
    if length(sacTime) > 0
        % downsample the pupil data
        pupilDN = downsample(pupil,downsampleRate);
        sacTimeInd = floor(sacTime/downsampleRate)+1;
        sacTimeInd = sacTimeInd(sacTimeInd>lookBeforeOnset) - lookBeforeOnset;
        
        % make a convolution matrix with width (putative IRF length)
        len = length(pupilDN);
        saccadeTimeVector = zeros(1, len);
        saccadeTimeVector(sacTimeInd) = 1;
        
        for trNum = 1:s.task{1}.trialnum-1
            sacDuringThisTr = sacTimeInd(sacTimeInd < edf.inds(trNum,2) & sacTimeInd > edf.inds(trNum,1));
        end
        
        d.putativeIRFtime = 4; % Important: how long you think IRF is.
        putativeIRFlength = d.putativeIRFtime*sampleRate/downsampleRate;
        
        for ii=1:putativeIRFlength
            Sacmatrix(1:len,ii) = [zeros(1,ii-1) saccadeTimeVector(1:len-ii+1)]';
        end
        
        Sacmatrix = Sacmatrix(1:len,:);
        pupilSacIrf = pupilDN' * pinv(Sacmatrix');
        sacIrf = [sacIrf; pupilSacIrf];
    end
    
    pupTrace{iRun} = pupil;
    
    %% estimate saccade rate function
    
    maxLength = max(edf.inds(:,2) - edf.inds(:,1)+1);
    rateRaster = nan(s.task{1}.trialnum, maxLength);
    saccadeTimeVector = zeros(1, length(pupil));
    saccadeTimeVector(sacTime) = 1;
    
    for trNum = 1:s.task{1}.trialnum-1
        rateRaster(trNum,1:edf.inds(trNum,2) - edf.inds(trNum,1)+1) = saccadeTimeVector(edf.inds(trNum,1):edf.inds(trNum,2));
    end
    
    if iRun>1
        if size(sacRate,2)<size(rateRaster,2)
            sacRate =resizeToMatch(sacRate,size(rateRaster,2),2);
        else
            rateRaster =resizeToMatch(rateRaster,size(sacRate,2),2);
        end
    end
    
    sacRate = [sacRate ; nanmean(rateRaster)];
    pupilAvg = [pupilAvg; nanmean(e.eye.pupil(:,1:floor(maxLength/100)*100))];
    baseline = [baseline; edf.gaze.baseline];

    keyboard
end

d.sacRate = sacRate;
d.sacIrf = sacIrf;
d.pupilAvg = pupilAvg + baseline; % add back baseline
d.pupil = pupTrace;
d.baseline = baseline;
d.sampleRate = sampleRate;
d.downsampleRate = downsampleRate;
d.sacTimes = sacTimes;
d.trInds = trInds;
d.accuracyByTrial = accuracyByTrial;

return;

