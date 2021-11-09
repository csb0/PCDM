function [samples] = randp (p, num)
    % Warrning! this function works when p is defined for 0-N, and each
    % index correspond to an integer (count of spikes). It is not designed 
    % to work for k-N for example.
    sampProbs = rand(num, 1); 
    samples = zeros(num, 1); % Pre-allocating
    for ii=1:length(p)-1
      criteria = sum(p(1:ii));
      samples(sampProbs>criteria) = ii; % It is definately not the efficient way, it overwrite on previous thing until it is right. However, it works.
                                        % At the end, it samples desired count of numbers and then in every correspondeing sample from                                        
                                        % the poisson distribution, it puts relative number from 0-20 based on cumulative poisson distribution
    end
end

