function [results, error] = fit2E(CpPredicted, times, singleBloodDraw, ...
    bloodDrawTime, singleBloodDrawErrFactor)
%FIT2E Summary of this function goes here
%   Detailed explanation goes here
    
    phi1 = @(b) trapz((CpPredicted - model2E(b,times)).^2) + singleBloodDrawErrFactor * ...
        ( (model2E(b, bloodDrawTime + 0.05) - model2E(b, bloodDrawTime - .05)) * 10 ...
        - singleBloodDraw) ^ 2;
    % Takes single blood draw into account for error.

    tic;
    results = simulatedAnnealing(phi1,...
        [0.1 0.002 0.1 0 0], ...
        [-1e6 0 -1e6 0 0], ...
        [0 0.01 0 0.001 2e6 ], 1e6, 1e-7, 0.90);
    toc;
    
    error = phi1(results);%/trapz(Cpint1.^2) ;
end

