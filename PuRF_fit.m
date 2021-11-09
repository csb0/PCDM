function [ RsquareFit ] = PuRF_fit(n,tMin,myIRF,sampleRate)
% PuRF_fit.m 
%       Authors: Charlie Burlingham, Saghar Mirbagheri
%
%       Date: Sept 26, 2020
%
%       Purpose: fits parametric form to saccade-locked pupil response
%
%       Inputs: - n, exponent that controls width of trough of Gamma-Erling function
%               - tMin, controls time-to-minimum of Gamma-Erling function
%               - myIRF, saccade-locked pupil response
%               - sampleRate, samping rate of saccade-locked pupil response

    timeStep = 1000/sampleRate;
    impulseTime = 4000;
    time = 0:timeStep:impulseTime-1;
    impulseForm = -1.*(time.^n).*exp(-n.*time./tMin);        
    MinIRF = min(myIRF);
    minSimIRF = min(impulseForm);
    ratioIRF = MinIRF /  minSimIRF;
    impulseForm =  ratioIRF * impulseForm;
    RsquareFit = mean((myIRF- (impulseForm) ).^2); 
end

