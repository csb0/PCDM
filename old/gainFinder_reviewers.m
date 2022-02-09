function [gain,genFunc] = gainFinder(saccRate,linearFilter,pupilTimeSeries,trialInds)
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
%            - saccRate: saccade rate function (1xN matrix, where 
%              N is the number of time samples within a trial)             
%            - linearFilter: estimated linear filter, based on fitting a
%              Gamma-Erlang function to the deconvolved saccade-locked
%              pupil response (1XN matrix)
%            - pupilTimeSeres: full pupil time series across all trials,
%              after blink interpolation and low-pass filtering
%              (1XCN matrix, where C is number of trials in a run)
%            - trialInds: indices of trial start and end times that
%              correspond to pupilTimeSeries (CX2 matrix)
%
%     Outputs: 
%            - gain: estimated gain 
%            - genFunc: estimate generator function (1XN matrix)

%% Estimate expected value of generator function for each trial type
threshold = 1; % threshold is arbitrary, but must be consistent
saccRate(saccRate==0) = eps; % if rate is zero set it to a very small number so we don't reach inifinity at end of gaussian tail
genFunc = 1-saccRate;
genFunc = norminv(genFunc,0,1);
genFunc  = threshold - genFunc;
genFunc = genFunc-nanmean(genFunc); % mean-subtract generator function

%% Compute time series of estimated generator function input across all trials (generatorSeries) 
trialLengthToFit = 2000; % 4 s at 500 Hz sampling rate, for example
generatorSeries = zeros(length(pupilTimeSeries),1);

for jj=1:size(trialInds,1)
    generatorSeries(trialInds(jj,1):trialInds(jj,2)) = genFunc(1:trialInds(jj,2)- trialInds(jj,1)+1);
end
generatorSeries = generatorSeries - mean(generatorSeries);

%% Predict pupil size by convolving generator series and linear filer
linearFilter = linearFilter-linearFilter(1); % linear filter should start from 0
pred = conv(linearFilter,generatorSeries);
predMatrix = pred(1:length(pupilTimeSeries));

%% Truncating trial lengths to desired lengths
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

%% Compute trial-averaged pupil size for data and model
trialPupilM = nan(size(trialInds,1), trialLengthToFit);
trialPupilD = nan(size(trialInds,1), trialLengthToFit);
for ii =1:size(trialInds,1)
    trialPupilM(ii,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = predMatrix(newTrInds(ii,1):newTrInds(ii,2))';
    trialPupilD(ii,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = pupilTimeSeries(newTrInds(ii,1):newTrInds(ii,2))';
end

predMatrixAvg = nanmean(trialPupilM);
predMatrixAvg = predMatrixAvg- (mean(predMatrixAvg,2)); % compute trial-averaged pupil size from model
pupilAvg = nanmean(trialPupilD); % compute trial-averaged pupil size from data

%% Solve for gain 
DM = [ones(length(pupilAvg),1), predMatrixAvg'];
sol = regress(pupilAvg',DM);
gain = sol(2:end);

end