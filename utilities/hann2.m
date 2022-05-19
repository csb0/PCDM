function [ window ] = hann2(windowLength)
%   Charlie Burlingham
%   May 19 2022
%
%   Purpose: creates a hanning window, works like hanning.m in Matlab's signal processing toolbox
%
%   Input:  windowLength, number of samples in the hanning window
%
%   Output: window, the hanning window

N = windowLength+2; % to match matlab hanning function (first and last samples of window removed
base = linspace(0,N,N);
window =  0.5*(1 - cos(2*pi*base/N));
window = window(2:end-1); % to match matlab hanning function (first and last samples of window removed
