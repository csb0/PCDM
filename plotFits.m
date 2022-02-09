function plotFits(d,f)
% plotFits.m
%
%     Cite: Burlingham C*, Mirbagheri S*, Heeger DJ (2022). Science 
%           Advances. *Equal Authors
%
%     Date: 2/9/21
%
%     Purpose: plots data vs. model fits
%

%% plot data vs. model fit

figure;
for ii = 1:4%length(f.gain)
    subplot(1,4,ii);
    tBase = linspace(0,(length(d.TEPR{ii})*d.downsampleRate)/d.sampleRate,length(d.TEPR{ii}));
    plot(tBase,d.TEPR{ii},'b','lineWidth',2);
    hold on; plot(tBase,f.pred{ii}(1:length(d.TEPR{ii}))+f.offset{ii},'k','lineWidth',2);
    box off; grid on;
    xlabel('Time (ms)')
    ylabel('Pupil area (AU)')
    title('Model fits')
end
legend('data','model')


end