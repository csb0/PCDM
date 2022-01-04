function [d] = dataAnalysisToolbox(in,op)

% todo - dec 28: make work with trials of different lengths (jittered
% timing) and with specifying the event to lock the pupil response to
% (e.g., button press, trial onset, etc.).




% usage: for use on a single run/block of trials

% outline:

% - BP filtering of pupil data
% - blink interpolation (optional)
% - saccade detection
% - saccade-locked IRF devonvolution
% - compute task-evoked pupil response  (user specifies event they want to lock pupil response to (task onset, or button press).
% - compute saccade rate function
% - compute offset

% inputs:

%   input struct with fields:
% xPos, yPos of eye (two vectors)
% pupil area (a vector)
% trial start time indexes (a vector)
% trial event times (a cell array where each cell is a block of trials, e.g. button press)
% sampling rate of eye tracker (Hz)
% trial types (e.g., easy / hard, or corr / incorr) - integers for each trial (e.g. 1,2,3,4,5)
% time window within a trial to make a prediction (e.g. 4 sec within a jittered ISI expt)

%   options struct with field:
% interpolateBlinks (0 or 1) to interpolate blinks from pupil time series (only use if you are inputting eyeLink data, otherwise, do this step yourself)
% downsampleRate for deconvolution of saccade-locked pupil response (default 5)

% outputs:

%

%% set options (move to another script eventually)
op.interpolateBlinks = 1;
op.downsampleRate = 5;

in.putativeIRFdur = 4; % how long you think IRF is, in seconds.

%% Interpolate blinks (Mathot's method, modified for Parker & Denison 2020)
if op.interpolateBlinks == 1
    in.pupilArea((isnan(in.pupilArea))) = 0; % set NaN regions to 0 for interp algorithm
    in.pupilArea2 = blinkinterp(in.pupilArea',in.sampleRate,5,4,50,75); % default options
end

%% Band-pass filter pupil signal
in.pupilArea3 = myBWfilter(in.pupilArea2,[0.03 , 10],in.sampleRate,'bandpass');

%% Detect small saccades and microsaccades (method of Engbert & Mergenthaler 2006)
vel = vecvel([in.xPos in.yPos], in.sampleRate, 3); % compute eye velocity
d.sacs = microsacc([in.xPos in.yPos], vel, 8, 7); % detect (micro)saccade times based on eye velocity

%% Estimate saccade-locked pupil response

% downsample the pupil data
sacOnset = d.sacs(:,1); % determine the indexes that saccades begin
pupilDN = downsample(in.pupilArea3,op.downsampleRate);
sacTimeInd = floor(sacOnset/op.downsampleRate)+1;

% make a convolution matrix with width (putative IRF length)
saccadeTimeVector = zeros(1, length(pupilDN)); saccadeTimeVector(sacTimeInd) = 1;
nTrials = length(in.startInds);
for trNum = 1:nTrials-1
    sacDuringThisTr = sacTimeInd(sacTimeInd < in.startInds(trNum,2) & sacTimeInd > in.startInds(trNum,1));
end
putativeIRFlength = in.putativeIRFdur*in.sampleRate/op.downsampleRate;
for ii=1:putativeIRFlength
    Sacmatrix(1:length(pupilDN),ii) = [zeros(1,ii-1) saccadeTimeVector(1:length(pupilDN)-ii+1)]';
end

Sacmatrix = Sacmatrix(1:length(pupilDN),:);
d.sacIRF = pupilDN * pinv(Sacmatrix'); % deconvolve

%% Estimate task-evoked pupil response  [add functionality to lock to whatever task event you want and for jittered timing]
nSamplesPerTrial = length(in.pupilArea)/nTrials;

in.pupilAreaMat = NaN(nTrials,nSamplesPerTrial);
for trNum = 1:nTrials-1
    curInds = in.startInds(trNum,1):in.startInds(trNum,2);
    in.pupilAreaMat(trNum,1:length(curInds)) = in.pupilArea3(curInds);
end

d.TEPR = nanmean(in.pupilAreaMat);

%% Estimate saccade rate function
rateRaster = nan(nTrials, nSamplesPerTrial);
saccadeTimeVector = zeros(1, length(in.pupilArea3));
saccadeTimeVector(sacOnset) = 1;

for trNum = 1:nTrials-1
    rateRaster(trNum,1:in.startInds(trNum,2) - in.startInds(trNum,1)+1) = saccadeTimeVector(in.startInds(trNum,1):in.startInds(trNum,2));
end

d.sacRate = nanmean(rateRaster);
d.sacRate = d.sacRate(1:nSamplesPerTrial);

%% Estimate additive offset
d.offset = nanmean(in.pupilArea);

end
