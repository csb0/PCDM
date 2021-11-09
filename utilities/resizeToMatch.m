% resizeToMatch
%   usage:resizeToMatch(signalToResize,idealSize,dimension)
%   by: Saghar
%   date: 06/28/2018
%   purpose: change the size of vector to match the length of a matrix so
%   it can get concatenated

function [resizedSignal] = resizeToMatch(signalToResize,idealSize,dimensionToWork)

if length(idealSize) == 1
    if dimensionToWork == 1
         if idealSize>size(signalToResize,1)
             resizedSignal = nan(idealSize,size(signalToResize,2));
             resizedSignal(1:size(signalToResize,1),:) = signalToResize;
        elseif idealSize<size(signalToResize,1)
             resizedSignal = signalToResize(1:idealSize,:); 
        else
             resizedSignal = signalToResize;     
        end
    elseif dimensionToWork == 2
        if idealSize>size(signalToResize,2)
             resizedSignal = nan(size(signalToResize,1),idealSize);
             resizedSignal(:,1:size(signalToResize,2)) = signalToResize;
        elseif idealSize<size(signalToResize,2)
             resizedSignal = signalToResize(:,1:idealSize); 
        else
             resizedSignal = signalToResize;     
        end
    else
        error('Please use 1 for working on the rows and 2 for columns');
    end
else
    error('Please use a number as size');
end

end