function [edf, newblinksmp] = myBlink_interpolate2(edf, plotme)
% interpolates blinks and missing data
% Anne Urai, 2016, modified by eli and charlie

% get the stuff we need
dat.time        = edf.gaze.time;
dat.pupil       = edf.gaze.pupil;
dat.gazex       = edf.gaze.x;
dat.gazey       = edf.gaze.y;
padding         = 0.150; % how long before and after do we want to pad?
data.fsample = edf.samplerate;
myNans = zeros(size(dat.pupil));
myNans2 = zeros(size(dat.pupil));
idx = [];
newPadblinksmp =[];
% for iBlink=2:length(edf.blinks.startTime)-2; % ATTN hack
for iBlink=1:length(edf.blinks.startTime) 
    if length(find(dat.time==edf.blinks.startTime(iBlink))) == 1  % if there is just one time that corresponds to the current blink. why is this here?
        if ~isempty(find(dat.time==edf.blinks.endTime(iBlink))) && iBlink~=length(edf.blinks.startTime) % CSB aug 12 2021: error handling, in case last blink end time is empty
            idx = cat(1, idx, [find(dat.time==edf.blinks.startTime(iBlink)) find(dat.time==edf.blinks.endTime(iBlink))]); % find indices of blink start and end times
        end
    end
end
blinksmp = idx;

% initialize settings
if ~exist('plotme', 'var'); plotme = true; end % plot all this stuff

% initialize output
newblinksmp = [];
if ~isempty(blinksmp),
    disp('interpolating EL-defined blinks');
    
    % ====================================================== %
    % STEP 1: INTERPOLATE EL-DEFINED BLINKS
    % ====================================================== %
    
    if plotme,
        clf;  sp1 = subplot(411); plot(dat.time,dat.pupil, 'color', [0.5 0.5 0.5]);
        axis tight; box off; ylabel('Raw');
        set(gca, 'xtick', []);
    end
    
    % merge 2 blinks into 1 if they are < 150 ms together (coalesce)
    coalesce = 0.150;
    for b = 1:size(blinksmp, 1)-1,
        if blinksmp(b+1, 1) - blinksmp(b, 2) < coalesce * data.fsample,
            blinksmp(b, 2) = blinksmp(b+1, 2);
            blinksmp(b+1, :) = nan;
        end
    end
    % remove those duplicates
    blinksmp(isnan(nanmean(blinksmp, 2)), :) = [];
    
    % pad the blinks
    padblinksmp(:,1) = round(blinksmp(:,1) - padding * data.fsample);
    padblinksmp(:,2) = round(blinksmp(:,2) + padding * data.fsample);
    
    % avoid idx outside range
    if any(padblinksmp(:) < 1), padblinksmp(find(padblinksmp < 1)) = 1; end
    if any(padblinksmp(:) > length(dat.pupil)), padblinksmp(find(padblinksmp > length(dat.pupil))) = length(dat.pupil); end
    
  %  keyboard
    % make the pupil NaN at those points
    for b = 1:size(padblinksmp,1),
        dat.pupil(padblinksmp(b,1):padblinksmp(b,2)) = NaN;
        dat.gazex(padblinksmp(b,1):padblinksmp(b,2)) = NaN;
        dat.gazey(padblinksmp(b,1):padblinksmp(b,2)) = NaN;
    end
    
    % also set the pupil to zero when there were missing data
    dat.pupil(dat.pupil < nanmedian(dat.pupil)-3*nanstd(dat.pupil)) = nan;
    dat.pupil(dat.pupil > nanmedian(dat.pupil)+3*nanstd(dat.pupil)) = nan;

    myNans = isnan(dat.pupil);
    
    % interpolate linearly
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');
    dat.gazex(isnan(dat.gazex)) = interp1(find(~isnan(dat.gazex)), ...
        dat.gazex(~isnan(dat.gazex)), find(isnan(dat.gazex)), 'linear');
    dat.gazey(isnan(dat.gazey)) = interp1(find(~isnan(dat.gazey)), ...
        dat.gazey(~isnan(dat.gazey)), find(isnan(dat.gazey)), 'linear'); 
  
    % to avoid edge artefacts at the beginning and end of file, pad in seconds
    % edgepad = 1;
    % dat.pupil(1:edgepad*data.fsample)           = NaN;
    % dat.pupil(end-edgepad*data.fsample : end)   = NaN;
    
    
    if plotme, sp2 = subplot(411); hold on;
        % show how well this worked
        plot(dat.time, dat.pupil, 'b');
        axis tight; box off; ylabel('Interp');
        set(gca, 'xtick', []);
    end
end
% also extrapolate ends
dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
    dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'nearest', 'extrap');
    
    

% ====================================================== %
% STEP 2: INTERPOLATE PEAK-DETECTED BLINKS
% ====================================================== %
% keyboard
% if any(isnan(dat.pupil))
%     return;
% end

assert(~any(isnan(dat.pupil)));
%win             = hanning(11);
%pupildatsmooth  = filter2(win.',dat.gazey,'same');

%dat.pupildiff   = (diff(pupildatsmooth) - mean(diff(pupildatsmooth))) / std(diff(pupildatsmooth));
%[peaks, loc]    = findpeaks(abs(dat.pupildiff), 'minpeakheight', 3*std(dat.pupildiff), 'minpeakdistance', 0.5*data.fsample);

dat.pupildiff = diff(dat.pupil);
aboveT = find(abs(dat.pupildiff) > 50);
coalesce = 0.150;

%peaks = [peaks1, peaks2];
%[loc, locOrder] = sort([loc1,loc2]);
%peaks = peaks(locOrder);


if plotme, sp2 = subplot(412);
    plot(dat.time(2:end), dat.pupildiff);
    hold on; plot(dat.time(loc), peaks, '.');
    box off; ylabel('Peaks');
    set(gca, 'xtick', []); ylim([-500 500]);
    sp3 = subplot(413); hold on;
end

if ~isempty(aboveT),
    newBlinksmp(:,1) = [aboveT(1), aboveT(find(diff(aboveT)>coalesce*data.fsample)+1)];
    newBlinksmp(:,2) = [aboveT(find(diff(aboveT)>coalesce*data.fsample)), aboveT(end)];



%keyboard
    if plotme,
        plot(dat.time, dat.pupil, 'color', [0.5 0.5 0.5]);
        hold on;
    end
    
    newPadblinksmp(:,1) = round(newBlinksmp(:,1) - padding * data.fsample);
    newPadblinksmp(:,2) = round(newBlinksmp(:,2) + padding * data.fsample);
    
    % avoid idx outside range
    if any(newPadblinksmp(:) < 1), newPadblinksmp(find(newPadblinksmp < 1)) = 1; end
    if any(newPadblinksmp(:) > length(dat.pupil)), newPadblinksmp(find(newPadblinksmp > length(dat.pupil))) = length(dat.pupil); end
    
    %keyboard
    for b = 1:size(newPadblinksmp,1),
        dat.pupil(newPadblinksmp(b,1):newPadblinksmp(b,2)) = NaN;
        dat.gazex(newPadblinksmp(b,1):newPadblinksmp(b,2)) = NaN;
        dat.gazey(newPadblinksmp(b,1):newPadblinksmp(b,2)) = NaN;
    end
    
    
    %keyboard
   % keyboard
    myNans2 = isnan(dat.pupil);
    % interpolate linearly
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');
     dat.gazex(isnan(dat.gazex)) = interp1(find(~isnan(dat.gazex)), ...
        dat.gazex(~isnan(dat.gazex)), find(isnan(dat.gazex)), 'linear');
     dat.gazey(isnan(dat.gazey)) = interp1(find(~isnan(dat.gazey)), ...
        dat.gazey(~isnan(dat.gazey)), find(isnan(dat.gazey)), 'linear');
end

%keyboard
% remove remaining nans (probably at the end)
dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
    dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'nearest', 'extrap');
newpupil = dat.pupil;

% link axes
if plotme,
    plot(dat.time, dat.pupil, 'b');
    ylim([min(dat.pupil)*0.9 max(dat.pupil)*1.1]);
    box off; ylabel('Clean');
    try
        linkaxes([sp1 sp2 sp3], 'x');
        set([sp1 sp2 sp3], 'tickdir', 'out');
    end
    xlim([dat.time(1)-10 dat.time(end)+10]);
end

allNans  = myNans+myNans2;
allNans = allNans>0;
%keyboard
% store the results
edf.gaze.nans = allNans;
edf.gaze.pupil = newpupil;
edf.gaze.x = dat.gazex;
edf.gaze.y = dat.gazey;
%edf.gaze.blinks = [padblinksmp; newblinksmp];
edf.gaze.blinks = [newPadblinksmp; padblinksmp];
end


