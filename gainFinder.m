function [d, f] = gainFinder(d,ii,op)
% gainFinder.m
%
%     Cite: Burlingham C*, Mirbagheri S*, Heeger DJ (2022). Science
%           Advances. *Equal Authors
%
%     Date: 2/9/22
%
%     Purpose: Solves for best-fit gain and generator function(s), and
%              returns model predictions and goodness of fit.
%

threshold = 1; % arbitrary, but fixed
d.sacRate2{ii}(d.sacRate2{ii}==0) = eps; % if rate is zero set it to a very small number so we don't reach inifinity at end of gaussian tail. Future implementations can instead regularize missing data..
Generator = 1-d.sacRate2{ii};
numTrialTypes = size(d.sacRate2{ii},1);
for kk = 1:numTrialTypes %% kk index is for trial types
    Generator(kk,:) = norminv(Generator(kk,:),0,1);
    Generator(kk,:)  = threshold - Generator(kk,:);
    Generator(kk,:) = Generator(kk,:)-nanmean(Generator(kk,:)); % mean-subtract generator function
end

MIR = d.parametricLinearFilter-d.parametricLinearFilter(1);

if op.fitTimeseries == 1
    
    % create generator series
    generatorSeries = zeros(length(d.pupilTS{ii}),max(d.trialTypes{ii}));
    predMatrix = zeros(length(d.pupilTS{ii}),max(d.trialTypes{ii}));
    for kk = 1:max(d.trialTypes{ii})
        thisType = find(d.trialTypes{ii}==kk);
        for jj=1:length(thisType) % looping through trials of that type
            generatorSeries(d.trInds{ii}(thisType(jj),1):d.trInds{ii}(thisType(jj),2),kk) = Generator(1,1:d.trInds{ii}(thisType(jj),2)- d.trInds{ii}(thisType(jj),1)+1)';
        end
        pred = conv(MIR,generatorSeries(:,kk));
        predMatrix(:,kk) = pred(1:length(d.pupilTS{ii}));
    end
    
    % remove bits beyond prediction window (e.g. 4 s past trial onset)
    preserveInds = [];
    for nn = 1:size(d.trInds{ii})
        if nn == size(d.trInds{ii},1)
            d.trInds{ii}(nn,2) = length(d.pupilTS{ii}); % make sure last trial index is equal to last index of pupil timeseries
            preserveInds = [preserveInds d.trInds{ii}(nn,1):d.trInds{ii}(nn,2)];
            preserveIndsMat(nn,:) = [d.trInds{ii}(nn,1) d.trInds{ii}(nn,2)];
        else
            preserveInds = [preserveInds d.trInds{ii}(nn,1):d.trInds{ii}(nn,1)+d.predictionWindow*d.sampleRate-1];
            preserveIndsMat(nn,:) = [d.trInds{ii}(nn,1) d.trInds{ii}(nn,1)+d.predictionWindow*d.sampleRate-1];
        end
    end
    predMatrix2 = predMatrix(preserveInds,:);
    d.pupilTS2{ii} = d.pupilTS{ii}(preserveInds);
    
    % fit
    DM = [ones(length(d.pupilTS2{ii}),1), predMatrix2];
    sol = regress(d.pupilTS2{ii}',DM);
    gain = sol(2:end);
    
    predTS = gain'*predMatrix2';
    f.predTS = predTS;
    
    % only works if all trials end up being cut to same length (i.e., won't
    % work in cases where some trials are shorter than prediction window -
    % update to below method **
    predAvgMat = NaN(size(d.trInds{ii},1),d.predictionWindow*d.sampleRate+10); % +10 is a buffer for fencepost issues
    for nn = 1:size(d.trInds{ii},1)
        startIndNew = find(preserveInds== preserveIndsMat(nn,1));
        endIndNew = find(preserveInds== preserveIndsMat(nn,2));
        predAvgMat(nn,1:endIndNew(end)-startIndNew(end)+1) = predTS(startIndNew(end):endIndNew(end));
    end
    pred = nanmean(predAvgMat(:,1:d.predictionWindow*d.sampleRate));
    pred = downsample(pred,d.downsampleRate,1);
    
    
    % cut trial-avg pupil data to length of prediction window
    d.TEPR{ii} = d.TEPR{ii}(1:d.predictionWindow*d.sampleRate/d.downsampleRate);
    
    offset = mean(d.TEPR{ii});
    
    SSres = sum((d.TEPR{ii}-(pred+offset)).^2);
    SStot = sum((d.TEPR{ii}-nanmean(d.TEPR{ii})).^2);
    Rsq = 1 - (SSres./SStot);
    
    
    %pred = nanmean(reshape(predTS',d.predictionWindow*d.sampleRate,length(d.trInds{ii}))');
    %pred = downsample(pred,d.downsampleRate,1);
    
    %{

trAvgPredMat = NaN(length(d.trInds{ii}),length(d.TEPR{ii}));
for ll = 1:length(d.trInds{ii})
    if ll == length(d.trInds{ii})
        trAvgPredMat(ll,1:length(pred(d.trInds{ii}(ll,1):end))) = pred(d.trInds{ii}(ll,1):end);
    else
        trAvgPredMat(ll,1:(d.trInds{ii}(ll,2)-d.trInds{ii}(ll,1))+1) = pred(d.trInds{ii}(ll,1): d.trInds{ii}(ll,2));
    end
end

pred = downsample(nanmean(trAvgPredMat),d.downsampleRate,1);
    %}
    
    
    
elseif op.fitTimeseries == 0
    
    for kk = 1:numTrialTypes     
        pupilAvg(kk,:) = d.TEPR_TT{ii}(kk,:);
        predT = cconv(MIR,Generator(kk,:),size(Generator,2));
        predT = predT- (mean(predT,2));
        predT2(kk,:) = downsample(predT,d.downsampleRate,1);
        
        DM = [ones(length(pupilAvg(kk,:)),1), predT2(kk,:)'];
        sol = regress(pupilAvg(kk,:)',DM);
        
        gain(kk) = sol(2);
        
        pred1(kk,:) = gain(kk).*predT2(kk,:); % final prediction
        pred(kk,:) = pred1(kk,1:d.predictionWindow*d.sampleRate/d.downsampleRate);
        
        d.TEPR_TT2{ii}(kk,:) = d.TEPR_TT{ii}(kk,1:d.predictionWindow*d.sampleRate/d.downsampleRate);
        
        offset(kk) = mean(d.TEPR_TT2{ii}(kk,:));
        
        SSres = sum((d.TEPR_TT2{ii}(kk,:)-(pred(kk,:)+offset(kk))).^2);
        SStot = sum((d.TEPR_TT2{ii}(kk,:)-nanmean(d.TEPR_TT2{ii}(kk,:))).^2);
        Rsq(kk) = real(1 - (SSres./SStot));
    end
end


% save out fits and parameters
f.gain = gain;
f.Generator = Generator;
f.offset = offset;
f.Rsq = Rsq;
f.pred = pred;

end

