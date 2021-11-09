function plotFits(d,f)
% plotFits.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 11/8/21
%
%     Purpose: plots data vs. model fits
%


%% plot data vs. model fit (trial-average responses)
figure;
for ii = 1:length(f.gain)
    subplot(1,length(f.gain),ii);
    offset = mean(d.pupilAvg(ii,:)) - mean(f.pred{ii});
    tBase = linspace(0,length(d.pupilAvg(ii,:))./d.sampleRate,length(d.pupilAvg(ii,:)));
    plot(tBase,d.pupilAvg(ii,:),'b','lineWidth',2);
    hold on; plot(tBase,f.pred{ii}(1:length(d.pupilAvg(ii,:)))+offset,'k','lineWidth',2);
    box off; grid on;
end
legend('data','model')


%% plot data vs. model fit (full time series)
% figure;
% ex = 5; % example run
%     tBase = linspace(0,length(f.predFullTS{ex})./d.sampleRate,length(f.predFullTS{ex}));
%     plot(tBase,f.predFullTS{ex})
%     hold on; plot(tBase,d.pupil{ex}(1:length(f.predFullTS{ex})))
%     box off; grid on;
% legend('data','model')

end