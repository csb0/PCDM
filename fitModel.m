function [d,f] = fitModel(directory)
% fitModel.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 11/8/21
%
%     Purpose: Fits Pupil Common Drive Model (PCDM) to pupil and saccade data
%
%     Usage:
%
%     Todo:
%              - build in conditional to check if edfToMat has already been run and if not run, run it

%% Estimate model inputs and parameters:
%  trial-average pupil response, linear filter, post-saccadic refractory period, saccade rate function, generator function, and gain

% go to directory with stimfiles and converted edf files
if ieNotDefined('directory'); directory = pwd; end
d = dataAnalysis(directory); % d is input structure

% estimate inter-saccadic interval and parameter k
[k, kCI] = fitK(d.sacTimes, d.sampleRate, 1);

% adjust saccade rate functions for estimated post-saccadic refractory period
d.sacRate2 = k.*d.sacRate;

% estimate parametric linear filter from run-average saccade-locked pupil response
IrfAvg = nanmean(d.sacIrf); % avg. saccade-locked IRF across all runs
[parametricLinearFilter, params, Rsq1, normFactor] = fitParametricPuRFfromData(IrfAvg,d.sampleRate./d.downsampleRate,d.sampleRate,1);
keyboard
% specify trial types (based on experimental design and/or analyses)
trialTypes = ones(1,length(d.trInds{1})); 


% find best gain 
for ii = 1:size(d.sacIrf,1) % loop across runs
    trialTypes = d.accuracyByTrial{ii}+1;
    d.pupil{ii} = d.pupil{ii}+d.baseline(ii);
    gain{ii} = gainFinder(nanmean(d.sacRate2),parametricLinearFilter,d.pupil{ii},trialTypes,d.trInds{ii});
end

% evaluate model with best gain
for ii = 1:size(d.sacIrf,1) % loop across runs
    [offset,pupilPrediction,trAvg,Rsq,Generator,predFullTimeSeries] = evalModel(gain(ii), nanmean(d.sacRate2),parametricLinearFilter,d.pupil{ii},trialTypes,d.trInds{ii});
    
    f.offset(ii) = offset; % additive offset
    f.gain{ii} = gain{ii}; % best-fit gain
    f.pred{ii} = pupilPrediction; % model prediction of task-evoked pupil response
    f.Rsq(ii) = Rsq; % R^2 of fit
    f.Generator(ii,:) = Generator; % mean-subtracted generator function
    f.avg(ii,:) = trAvg(1:2000); % trial-average task-evoked pupil response    % ATTN:harcoded trial length
    f.predFullTS{ii} = predFullTimeSeries;
end

%% plot data and model fits
plotFits(d,f);

keyboard

end