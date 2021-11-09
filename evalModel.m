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

threshold =  1;
saccRate(saccRate==0) = eps; % set zero rates to a very small number so we don't reach inifinity at end of gaussian tail
Generator = 1-saccRate;
Generator = norminv(Generator,0,1);
Generator  = threshold - Generator;
meanGen = nanmean(Generator);
Generator = Generator-meanGen; % mean-subtract generator function.

saccIrf = saccIrf-saccIrf(1); 

generatorSeries = zeros(length(pupil),max(trialType)); 
predMatrix = zeros(length(pupil),max(trialType)); 

trialLengthToFit = 2000;

for ii=1:max(trialType)
    thisType = find(trialType==ii);
    for jj=1:length(thisType)
        generatorSeries(trialInds(thisType(jj),1):trialInds(thisType(jj),2),ii) = Generator(ii,1:trialInds(thisType(jj),2)- trialInds(thisType(jj),1)+1);
    end
    pred = conv(saccIrf,generatorSeries(:,ii));
    predMatrix(:,ii) = pred(1:length(pupil));
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

pupil(inds) = [];
predMatrix(inds,:) =[];

trialPupilM = nan(size(trialInds,1),max(trialType), trialLengthToFit);
trialPupilD = nan(size(trialInds,1), trialLengthToFit);
for ii =1:size(trialInds,1)
     trialPupilM(ii,:,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = predMatrix(newTrInds(ii,1):newTrInds(ii,2),:);
     trialPupilD(ii,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = pupil(newTrInds(ii,1):newTrInds(ii,2),:);
end

predMatrixAvg = squeeze(nanmean(trialPupilM)); 
dataPupilAvg = nanmean(trialPupilD);

pred = predMatrixAvg*gain;
predFullTimeSeries = predMatrix*gain; % CSB


offset = (mean(dataPupilAvg)-mean(pred));
predPupilAvg = pred + offset;

offset = (mean(pupil)-mean(predFullTimeSeries)); % CSB
predFullTimeSeries = predFullTimeSeries + offset; %CSB

SSres = sum((dataPupilAvg'-predPupilAvg).^2);
SStot = sum((dataPupilAvg-mean(dataPupilAvg)).^2);
Rsq = 1 - (SSres./SStot);


end

