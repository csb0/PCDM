function [filteredSig] = myBWfilter(signal,cuttOff,sRate,type)

if strcmp(type,'bandpass') || strcmp(type,'stop')
    [b, a] = butter(1,cuttOff/(sRate/2),type);
else
    [b, a] = butter(2,cuttOff/(sRate/2),type);
end

signal = [repmat(signal(1,1),1,length(signal)), signal, flip(signal)];

filteredSignal =  filtfilt(b,a,signal);
filteredSig = filteredSignal(1+length(filteredSignal)/3:length(filteredSignal)/3*2);

end

