function plotFits(d,f)
% plotFits.m
%
%     Authors: Charlie S. Burlingham & Saghar Mirbagheri
%
%     Date: 2/8/21
%
%     Purpose: plots data vs. model fits
%
%     Todos: - improve this function, by having it plot everything for
%              diagnostic checks. 
%            - have it split by trial type             


%% plot data vs. model fit (trial-average responses)
figure;
for ii = 1:length(f.gain)
    subplot(1,length(f.gain),ii);
    tBase = linspace(0,length(d.TEPR{ii})*d.downsampleRate,length(d.TEPR{ii}));
    plot(tBase,d.TEPR{ii},'b','lineWidth',2);
    hold on; plot(tBase,f.pred{ii}(1:length(d.TEPR{ii}))+f.offset{ii},'k','lineWidth',2);
    box off; grid on;
    xlabel('Time (ms)')
    ylabel('Pupil area (AU)')
    title('Model fits')
end
legend('data','model')



end