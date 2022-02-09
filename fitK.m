function [k kCI] = fitK(saccs, samplingRate, plotFit)
% fitK.m
%
%     Cite: Burlingham C*, Mirbagheri S*, Heeger DJ (2022). Science 
%           Advances. *Equal Authors
%
%     Date: 2/9/21
%
%     Purpose: Fits k, the shape parameter of a gamma distribution, from  
%              the measured inter-saccadic interval statistics.
%
%              In our model, k controls the refractory period following 
%              each saccade by removing every k-th saccade (equivalent to
%              scaling the  trial-averaged rate function).   
%                  
%
%              Inputs: 
%               - saccs, a cell array containing one matrix with saccade 
%               onset, offset, amplitde, etc... for each run of trials.
%               Each matrix is what is returned by default by the Engbert &
%               Kliegel algorithm. trMSAnalMicLeft_Vs_Right.m or
%               dataAnalysis.m will output this cell array.
%               - samplingRate, sampling rate of eye tracker or input data.
%               500 Hz by default
%
%               Outputs:
%               - k, remove every k-th saccade to model the empirical 
%               inter-saccadic interval
%               - kCI, 95% confidence interval around estimate of k.

if ieNotDefined('samplingRate')
    samplingRate = 500; %Hz
end

if ieNotDefined('plotFit')
    plotFit = 1; % plot fit
end

% concatenate data from different runs
interSacIntervals = [];
for ii = 1:length(saccs)
    if ~isempty(saccs{ii})
        
        temp = saccs{ii}(:,1:2);
        temp(:,2) = circshift(temp(:,2),1);
        temp = [temp(2:end,2) temp(2:end,1) ]; % time difference between end of last saccade and beginning of current saccade
        diffTempSec = (diff(temp,[],2)./samplingRate).*1000; % ms
        
        % set a minimum interval based on what is biomechanically possible.
        % following Otero-Millan 2008, we used 20 ms
        diffTempSec(diffTempSec<20) = []; 
        
        
        interSacIntervals = [interSacIntervals; diffTempSec ]; % ms. compute intersaccade interval as finite difference of saccade onset times
    end
end

% fit a gamma distribution to the measured inter-saccade intervals
[phat pci] = gamfit(interSacIntervals);

k = phat(1);
theta = phat(2);

kCI = pci(:,1);

if plotFit
    % plot the fitted gamma distribution over the normalized histogram
    figure(1); subplot(2,1,1); histogram(interSacIntervals,400,'Normalization','pdf')
    hold on;
    plot(0:.01:max(interSacIntervals),gampdf(0:.01:max(interSacIntervals),k,theta),'r','lineWidth',2);
end

title('Diagnostic check: fit of Gamma distribution to inter-saccade interval distribution')

end