function plotFits(d,f,op)
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

if op.fitTimeseries == 0
    
    figure;
    counter = 0;
    for jj = 1:max(d.trialTypes{1})
        for ii = 1:length(f.gain)
            counter = counter+1;
            subplot(max(d.trialTypes{1}),length(f.gain),counter);
            tBase = linspace(0,(length(d.TEPR_TT2{ii}(jj,:))*d.downsampleRate)/d.sampleRate,length(d.TEPR_TT2{ii}(jj,:)));
            plot(tBase,d.TEPR_TT2{ii}(jj,:),'b','lineWidth',2);
            hold on; plot(tBase,f.pred{ii}(jj,:)+f.offset{ii}(jj),'k','lineWidth',2);
            box off; grid on;
            xlabel('Time (ms)')
            ylabel('Pupil area (AU)')
            title(['Run #' num2str(ii) ', trial type #' num2str(jj)])
        end
    end
    legend('data','model')
    
elseif op.fitTimeseries == 1
    
    
    for ii = 1:length(f.gain)
        figure(5)
        subplot(1,length(f.gain),ii);
        tBase = linspace(0,(length(d.TEPR{ii})*d.downsampleRate)/d.sampleRate,length(d.TEPR{ii}));
        plot(tBase,d.TEPR{ii},'b','lineWidth',2);
        hold on; plot(tBase,f.pred{ii}(1:length(d.TEPR{ii}))+f.offset{ii},'k','lineWidth',2);
        box off; grid on;
        xlabel('Time (ms)')
        ylabel('Pupil area (AU)')
        title(['Run #' num2str(ii)])
        
        figure(6)
        subplot(length(f.gain),1,ii);
        tBase = linspace(0,length(d.pupilTS2{ii})/d.sampleRate,length(d.pupilTS2{ii}));
        plot(tBase,d.pupilTS2{ii},'b','lineWidth',2);
        hold on; plot(tBase,f.predTS{ii},'k','lineWidth',2);
        box off; grid on;
        xlabel('Time (ms)')
        ylabel('Pupil area (AU)')
        title(['Run #' num2str(ii)])
        
    end
    legend('data','model')
    
end

end