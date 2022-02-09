function f = gainFinder(d,ii)
% gainFinder.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 11/8/21
%
%     Purpose:  Solves for best-fit gain and generator function(s), and
%               returns model predictions and goodness of fit.
%

threshold = 1; % arbitrary, but fixed 
d.sacRate2{ii}(d.sacRate2{ii}==0) = eps; % if rate is zero set it to a very small number so we don't reach inifinity at end of gaussian tail. Future implementations can instead regularize missing data..
Generator = 1-d.sacRate2{ii};
Generator = norminv(Generator,0,1);
Generator  = threshold - Generator;
Generator = Generator-nanmean(Generator); % mean-subtract generator function

MIR = d.parametricLinearFilter-d.parametricLinearFilter(1);

for ii = 1:size(Generator,1)
    predT = cconv(MIR,Generator(ii,:),size(Generator,2));
    predT = predT- (mean(predT,2));
    predT = downsample(predT,d.downsampleRate,1);
    predEach(ii,:) = predT;
    
    DM = [ones(length(d.TEPR{ii}),1), predEach(ii,:)'];
    sol = regress(d.TEPR{ii}',DM);
    
    gain(ii) = sol(2:end);
end

pred = gain*predEach; % final prediction

offset = (mean(d.TEPR{ii})-mean(pred));

SSres = sum((d.TEPR{ii}-pred).^2);
SStot = sum((d.TEPR{ii}-mean(d.TEPR{ii})).^2);
Rsq = 1 - (SSres./SStot);

f.gain = gain;
f.Generator = Generator;
f.offset = offset;
f.Rsq = Rsq;
f.pred = pred;

end

