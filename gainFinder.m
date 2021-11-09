function [gain] = gainFinder(myRate,MIR,pupil,trialType,trialInds)
% gainFinder.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 11/8/21
%
%     Purpose:  This can be used with fminsearch to fit the best-fit gain.
%
%     Usage:
%

threshold = 1; % arbitrary
myRate(myRate==0) = eps; % if rate is zero set it to a very small number so we don't reach inifinity at end of gaussian tail
Generator = 1-myRate;
Generator = norminv(Generator,0,1);
Generator  = threshold - Generator;
Generator = Generator-nanmean(Generator); % mean-subtract generator function

trialLengthToFit = 2000;
MIR = MIR-MIR(1); 
generatorSeries = zeros(length(pupil),max(trialType)); 
predMatrix = zeros(length(pupil),max(trialType)); 
for ii=1:max(trialType)
    thisType = find(trialType==ii);
    for jj=1:length(thisType)
        % CSB: need to estimate multiple generator function here, one per
        % trial type ? or use just one genFunc for all trialTypes
        generatorSeries(trialInds(thisType(jj),1):trialInds(thisType(jj),2),ii) = Generator(1,1:trialInds(thisType(jj),2)- trialInds(thisType(jj),1)+1);
    end
    generatorSeries = generatorSeries - mean(generatorSeries(:,ii));
    pred = conv(MIR,generatorSeries(:,ii));
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
     trialPupilM(ii,:,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = predMatrix(newTrInds(ii,1):newTrInds(ii,2),:)';
     trialPupilD(ii,1:newTrInds(ii,2)-newTrInds(ii,1)+1) = pupil(newTrInds(ii,1):newTrInds(ii,2),:)';
end
keyboard
predMatrixAvg = squeeze(nanmean(trialPupilM)); 
predMatrixAvg = predMatrixAvg- (mean(predMatrixAvg));
pupilAvg = nanmean(trialPupilD);

DM = [ones(length(pupilAvg),1), predMatrixAvg'];
sol = regress(pupilAvg',DM);

gain = sol(2:end);

end

