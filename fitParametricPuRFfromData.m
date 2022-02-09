function [parametricLinearFilter, params, Rsq, normFactor] = fitParametricPuRFfromData(irfData,samplingRate,outSamplingRate,plotFit)
% fitParametricPuRFfromData.m 
%
%     Cite: Burlingham C*, Mirbagheri S*, Heeger DJ (2022). Science 
%           Advances. *Equal Authors
%
%     Date: 2/9/22
%
%     Purpose: Fits a parametric form to the observed saccade-locked pupillary
%             impulse response function. This is a poor alternative to un-deforming the
%             IRF according to our model, but is better than using a fixed parametric
%             form for the PuRF, which is what everyone has done previously. This at
%             least captures the inter-individual variation in the time-to-minimum and
%             the the width of the IRF around its minimum point.
%
%     Inputs: - irfData, a vector containing the saccade-locked pupil response
%             - samplingRate, input sampling rate of the data, either from eyeTracker directly or after downsampling
%             - plotFit, 0 or 1, where 1 plots the measured and fit IRFs on top of eachother
%             - outSamplingRate, the sampling rate you want the parametric filter to have
%
%     Outputs - parametricLinearFilter, the parametric fit to the PuRF, a gamma-erling function
%             - params, the params of the Erling function (n and tMin)
%             - Rsq, the R^2 of the fit

%% Fit parametric form
irfData = irfData-irfData(1);
objFunc = @(n,tMin) PuRF_fit(n,tMin,irfData,samplingRate); % returns rms error between measured IRF and gamma-erlang parametric form

x0 = [rand()*20 (2000-300)*rand()+300]; % randomize initial points for params

params = fminsearch(@(x) objFunc(x(1),x(2)),x0); % optimize.
% params(1) is n
% params(2) is tMin
%params(1) = params(1).*2; % uncomment if you want to adjust width to "undeform" in a VERY adhoc way
timeStep = 1000/outSamplingRate;
impulseTime = 4000; % 4 seconds. this is an assumption used throughout, based on observation that saccade-locked purf goes back to baseline after 4 sec
time = 0:timeStep:impulseTime-1;
parametricLinearFilter = -1.*(time.^params(1)).*exp(-params(1).*time./params(2)); % parametric form. gamma-erlang

% just for compuing R^2, also compute the filter in the input timebase
timeStep = 1000/samplingRate;
impulseTime = 4000; % 4 seconds. this is an assumption used throughout, based on observation that saccade-locked purf goes back to baseline after 4 sec
time = 0:timeStep:impulseTime-1;
parametricLinearFilterOrginalSamplingRate = -1.*(time.^params(1)).*exp(-params(1).*time./params(2)); % parametric form. gamma-erlang


%% Normalize the PuRF amplitude just for plotting. we will later renormalize it when we track gain over time
normFactor = abs(min(parametricLinearFilter))./abs(min(irfData));
parametricLinearFilter = parametricLinearFilter./normFactor; % nov 11, 2020: checked and it barely affects best-fit sigma value (neglibly) to normalize this or not

% compute R^2 of fit
parametricLinearFilterOrginalSamplingRate = parametricLinearFilterOrginalSamplingRate./normFactor;
SSres = sum((irfData-parametricLinearFilterOrginalSamplingRate).^2);
SStot = sum((irfData-mean(irfData)).^2);
Rsq = 1 - (SSres./SStot);

%% plot fit vs data optionally

if plotFit
     tBase = linspace(0,impulseTime/1000,length(parametricLinearFilterOrginalSamplingRate));
     figure(1); subplot(2,1,2);
     hold on;
     plot(tBase,irfData,'k','lineWidth',2);
     plot(tBase,parametricLinearFilterOrginalSamplingRate,'r','lineWidth',2);
     xlabel('Time (s)')
     ylabel('Pupil size (AU)')
     legend('data','parametric fit')
     box off; grid on;
end
title('Diagnostic check: fit of linear filter to saccade-locked pupil response')

end
