function [d] = dataAnalysis(in)
% dataAnalysis.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 2/8/21
%
%     Purpose: Pre-processes eye data and estimates inputs to model.
%
%              Specifically does the following: blink interpolation
%              (optional), bandpass filtering of pupil data, saccade
%              detection, devonvolution of saccade-locked IRF, estimation
%              of task-evoked pupil response (user can format input data
%              to lock pupil response to various events, e.g., task onset,
%              cues, or button press), estimation of saccade rate
%              function(s) (one per trial type, e.g., correct vs. error
%              trials).
%
%     Inputs:
%
%              input struct "in" with fields containing cell array with one
%              cell per run of data. Field names:
%
%              - yPos: horizontal gaze position (column vector)
%              - yPos: vertical gaze position (column vector)
%              - pupilArea: pupil area (column vector)
%              - startInds: trial start indexes in samples (n x 2 matrix
%                with trial start and end times, n is number of trials)
%              - sampleRate: sampling rate of eye tracker (Hz)
%              - trialTypes: trial types (e.g., easy / hard, or corr\error)
%                integers for each trial, e.g. 1,2,3,4,5 (1 x n vector)
%              - predictionWindow: time window within a trial to make a
%                prediction (e.g. 4 sec within a jittered ISI expt)
%
%              options struct "op" with field names:
%
%              - interpolateBlinks (0 or 1) to interpolate blinks from pupil
%                time series (only use if you are inputting eyeLink data,
%                otherwise, do this step yourself)
%              - downsampleRate for deconvolution of saccade-locked pupil
%                response (default 5)
%

%% Set options
op.interpolateBlinks = 1; % 1, blink interpolation on; 0, interpolation off (use this if your input data is already blink interpolated)
op.downsampleRate = 5; % downsample rate for the deconvolution design matrix. Higher numbers makes code run faster, but shouldn't be set higher than Nyquist frequency.

in.putativeIRFdur = 4; % how long you think the saccade-locked pupil response is in seconds. User-defined, but 4 s is a good estimate according to our results and the literature
%in.trialTypes = ones(1,size(in.startInds,1)); % by default set to all ones. Change this to fit different trial types. For example, set all correct trials to 1 and all error trials to 2 if you want to estimate arousal on correct vs. incorrect trials.
% debug: to generate random trial types
%in.trialTypes( randi(75,75,1)) = 2; % CSB: remove when not using ***

in.predictionWindow = 4; % The time window you want to make the prediction for the TEPR. Usually the max trial length, but can be shorter if you want.


%%
numRuns = length(in.pupilArea); % number of runs of data (number of cells passed as input)

for kk = 1:numRuns
    
    %% Interpolate blinks (Mathot's method, modified for Parker & Denison 2020)
    if op.interpolateBlinks == 1
        in.pupilArea{kk}((isnan(in.pupilArea{kk}))) = 0; % set NaN regions to 0 for interp algorithm
        in.pupilArea2 = blinkinterp(in.pupilArea{kk}',in.sampleRate{kk},5,4,50,75); % default options
    end
    
    %% Band-pass filter pupil signal
    in.pupilArea3 = myBWfilter(in.pupilArea2,[0.03 , 10],in.sampleRate{kk},'bandpass');
    
    %% Detect small saccades and microsaccades (method of Engbert & Mergenthaler 2006)
    vel = vecvel([in.xPos{kk} in.yPos{kk}], in.sampleRate{kk}, 3); % compute eye velocity
    d.saccTimes{kk} = microsacc([in.xPos{kk} in.yPos{kk}], vel, 8, 7); % detect (micro)saccade times based on eye velocity
    
    %% Estimate saccade-locked pupil response
    
    % downsample the pupil data
    sacOnset = d.saccTimes{kk}(:,1); % determine the indexes that saccades begin
    pupilDN = downsample(in.pupilArea3,op.downsampleRate);
    sacTimeInd = floor(sacOnset/op.downsampleRate)+1;
    
    % make a convolution matrix with width (putative IRF length)
    saccadeTimeVector = zeros(1, length(pupilDN)); saccadeTimeVector(sacTimeInd) = 1;
    nTrials = length(in.startInds);
    for trNum = 1:nTrials-1
        sacDuringThisTr = sacTimeInd(sacTimeInd < in.startInds{kk}(trNum,2) & sacTimeInd > in.startInds{kk}(trNum,1));
    end
    putativeIRFlength = in.putativeIRFdur*in.sampleRate{kk}/op.downsampleRate;
    for ii=1:putativeIRFlength
        Sacmatrix(1:length(pupilDN),ii) = [zeros(1,ii-1) saccadeTimeVector(1:length(pupilDN)-ii+1)]';
    end
    
    Sacmatrix = Sacmatrix(1:length(pupilDN),:);
    d.sacIrf(kk,:) = pupilDN * pinv(Sacmatrix'); % deconvolve
    
    %% Estimate task-evoked pupil response
    
    maxTrialLength = max(in.startInds{kk}(:,2)-in.startInds{kk}(:,1))+1; % samples
    
    % make a convolution matrix with width (max trial length)
    trialStartTimeDN = floor(in.startInds{kk}(:,1)/op.downsampleRate)+1;
    trialStartTimeVec = zeros(1, length(pupilDN)); trialStartTimeVec(trialStartTimeDN) = 1;
    
    maxTrialLengthDN = maxTrialLength/op.downsampleRate;
    for ii=1:maxTrialLengthDN
        TEPRmatrix(1:length(pupilDN),ii) = [zeros(1,ii-1) trialStartTimeVec(1:length(pupilDN)-ii+1)]';
    end
    
    TEPRmatrix = TEPRmatrix(1:length(pupilDN),:);
    d.TEPR{kk} = pupilDN * pinv(TEPRmatrix'); % deconvolve
    
    %% Estimate saccade rate function(s) (one per trial type)
    
    in.numTrialTypes{kk} = max(in.trialTypes{kk}); % number of trial types within the run of data. Shouldn't be set too high without feeding in much more data to constrain estimation.
    
    for ii = 1:double(in.numTrialTypes{kk}) % loop over trial types
        trialNumsPerType{ii} = find(in.trialTypes{kk}(1:end-1)==ii);
        nTrialsPerType(ii) = length(trialNumsPerType{ii});
        rateRaster{ii} = nan(nTrialsPerType(ii), maxTrialLength);
        saccadeTimeVector2{ii} = zeros(1, length(in.pupilArea3));
        saccadeTimeVector2{ii}(sacOnset) = 1;
        
        for trNum = trialNumsPerType{ii}
            rateRaster{ii}(trNum,1:in.startInds{kk}(trNum,2) - in.startInds{kk}(trNum,1)+1) = saccadeTimeVector2{ii}(in.startInds{kk}(trNum,1):in.startInds{kk}(trNum,2));
        end
        
        d.sacRateTemp = nanmean(rateRaster{ii});
        d.sacRate{kk}(ii,:) = d.sacRateTemp(1:maxTrialLength);
    end
    
end

%% Save out other variables
d.trialTypes = in.trialTypes;
d.sampleRate = in.sampleRate{1}; % there should be only one sampling rate for all data, or remove the {1}
d.downsampleRate = op.downsampleRate;
d.trInds = in.startInds;
d.predictionWindow = in.predictionWindow;



end