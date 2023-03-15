function [d,f] = fitModel(in)
% fitModel.m
%
%     Cite: Burlingham C*, Mirbagheri S*, Heeger DJ (2022). Science
%           Advances. *Equal Authors
%
%     Date: 2/9/21
%
%     Purpose: Fits Pupil Common Drive Model (PCDM) to pupil and saccade
%              data. Estimates model inputs and parameters: trial-average
%              pupil response, linear filter, post-saccadic refractory
%              period, saccade rate function, generator function, and gain.
%
%     Usage:   Input your eye data and task/contrast info (see
%              dataAnalysis.m for format and details).
%
%%

% Load in data struct "in" (see dataAnalysis for format and details)
[d, op] = dataAnalysis(in); % d is input structure

% estimate inter-saccadic interval and parameter k
[k, k_CI] = fitK(d.saccTimes, d.sampleRate, 1);

% adjust saccade rate functions for estimated post-saccadic refractory period
numRuns = size(d.sacIrf,1); % loop across runs
for ii = 1:numRuns
    d.sacRate2{ii} = round(k).*d.sacRate{ii};
end

% estimate parametric linear filter from run-average saccade-locked pupil response
IrfAvg = nanmean(d.sacIrf,1); % avg. saccade-locked IRF across all runs
[parametricLinearFilter, params_PuRF, Rsq_PuRF, normFactor_PuRF] = fitParametricPuRFfromData(IrfAvg,d.sampleRate./d.downsampleRate,d.sampleRate,1);
d.parametricLinearFilter = parametricLinearFilter;

% fit and evaluate model parameters
for ii = 1:numRuns % loop across runs
    [d, temp] = gainFinder(d,ii,op); % ii index is for runs
    
    f.gain{ii} = temp.gain; % gain
    f.Generator{ii} = temp.Generator; % generator function
    f.offset{ii} = temp.offset; % additive offset
    f.Rsq{ii} = temp.Rsq; % R-squared of model
    f.pred{ii} = temp.pred; % model prediction of pupil response
    
    if op.fitTimeseries == 1
        f.predTS{ii} = temp.predTS; % model prediction of full timeseries
    end
end

% plot data and model fits
plotFits(d,f,op);

end