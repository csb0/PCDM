function [gain,genFunc] = gainFinder(saccRate,linearFilter,pupilTimeSeries,trialType,trialInds)
% gainFinder.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date:    11/15/21
%
%     Purpose: estimates gain and generator function from pupil size,
%              trial-averaged saccade rate, and estimated linear filter 
%
%     Inputs: 
%            - saccRate: saccade rate function (MxN matrix, where M is
%              number of trial types, and N is the number of time samples
%              within a trial)
%            - linearFilter: estimated linear filter, based on fitting a
%              Gamma-Erlang function to the deconvolved saccade-locked
%              pupil response (1XN matrix)
%            - pupilTimeSeres: full pupil time series across all trials,
%              after blink interpolation and low-pass filtering
%              (1XCN matrix, where C is number of trials in a run)
%            - trialType: trial type for each trial in a run (e.g., 1 or 2 
%              for correct or incorrect trials). Works for any contrast 
%              you want (1XC matrix)
%            - trialInds: indices of trial start and end times that
%              correspond to pupilTimeSeries (CX2 matrix)
%
%     Outputs: 
%            - gain: estimated gain for each trial type 
%              (1XM matrix)
%            - genFunc: estimated generator function for each trial type
%              (MXN matrix)
%

%% Estimate expected value of generator function for each trial type
threshold = 1; % threshold is arbitrary, but must be consistent
for ii = 1:size(saccRate,1)
    curSaccRate = saccRate(ii,:);
    curSaccRate(curSaccRate==0) = eps; % if rate is zero, set it to a very small number epsilon so we don't reach inifinity at end of gaussian tail
    curGenFunc = 1-curSaccRate;
    curGenFunc = norminv(curGenFunc,0,1);
    curGenFunc  = threshold - curGenFunc;
    curGenFunc = curGenFunc-nanmean(curGenFunc); % mean-subtract generator function
    genFunc(ii,:) = curGenFunc;
end

%% Compute time series of estimated generator function input across all trials (generatorSeries) 
trialLengthToFit = 2000; % 4 s at 500 Hz sampling rate, for example
linearFilter = linearFilter-linearFilter(1);
generatorSeries = zeros(length(pupilTimeSeries),max(trialType));
predMatrix = zeros(length(pupilTimeSeries),max(trialType));
for ii=1:max(trialType)
    thisType = find(trialType==ii);
    for jj=1:length(thisType) 
        generatorSeries(trialInds(thisType(jj),1):trialInds(thisType(jj),2),ii) = genFunc(ii,1:trialInds(thisType(jj),2)- trialInds(thisType(jj),1)+1);
    end
end
generatorSeries = generatorSeries - mean(mean(generatorSeries));

%% Predict pupil size by convolving generator series and linear filer
for ii=1:max(trialType)
    pred = conv(linearFilter,generatorSeries(:,ii));
    predMatrix(:,ii) = pred(1:length(pupilTimeSeries));
end

inds =[];
newTrInds = trialInds;
for ii=1:size(trialInds,1)
    if (trialInds(ii,2)-trialInds(ii,1)) >  trialLengthToFit
        inds = [inds, trialInds(ii,1)+ trialLengthToFit:trialInds(ii,2)];
        newTrInds(ii,2) =  newTrInds(ii,1)+ trialLengthToFit - 1;
        newTrInds(ii+1:end,:) =  newTrInds(ii+1:end,:) - length(trialInds(ii,1)+ trialLengthToFit:trialInds(ii,2));
    end
end
pupilTimeSeries(inds) = [];
predMatrix(inds,:) =[];

trialPupilM = nan(size(trialInds,1),max(trialType), trialLengthToFit);
trialPupilD = nan(size(trialInds,1), trialLengthToFit);
for ii =1:size(trialInds,1)
    trialPupilM(ii,:,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = predMatrix(newTrInds(ii,1):newTrInds(ii,2),:)';
    trialPupilD(ii,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = pupilTimeSeries(newTrInds(ii,1):newTrInds(ii,2),:)';
end

%% Compute trial-averaged pupil size for data and model
predMatrixAvg = squeeze(nanmean(trialPupilM));
predMatrixAvg = predMatrixAvg- (mean(predMatrixAvg,2)); % compute trial-averaged pupil size from model
pupilAvg = nanmean(trialPupilD); % compute trial-averaged pupil size from data

%% Solve for gain 
DM = [ones(length(pupilAvg),1), predMatrixAvg'];
sol = regress(pupilAvg',DM);
gain = sol(2:end);
keyboard
end