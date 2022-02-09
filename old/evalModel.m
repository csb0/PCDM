function [offset,predPupilAvg,dataPupilAvg,Rsq,Generator,predFullTimeSeries] = evalModel(gain,saccRate,saccIrf,pupil,trialType,trialInds)
% evalModel.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 11/8/21
%
%     Purpose: This evaluates model with the best fit sigma and extracts the
%          offset and trial-average pupil prediction.
%
%     Usage:
%

threshold = 1; % arbitrary
d.sacRate2(d.sacRate2==0) = eps; % if rate is zero set it to a very small number so we don't reach inifinity at end of gaussian tail
Generator = 1-d.sacRate2;
Generator = norminv(Generator,0,1);
Generator  = threshold - Generator;
Generator = Generator-nanmean(Generator); % mean-subtract generator function

MIR = d.parametricLinearFilter-d.parametricLinearFilter(1);



pred = predMatrixAvg*gain;


offset = (mean(dataPupilAvg)-mean(pred));
predPupilAvg = pred + offset;


SSres = sum((dataPupilAvg'-predPupilAvg).^2);
SStot = sum((dataPupilAvg-mean(dataPupilAvg)).^2);
Rsq = 1 - (SSres./SStot);


end

