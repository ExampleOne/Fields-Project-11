function [ optSoFar ] = simulatedAnnealing( objFunc, initial, lower, ...
    upper, temperature, threshold, coolingRate)
% SIMULATEDANNEALING Optimises an objective function using simulated
% annealing.
%   Uses simulated annealing to optimise objective function objFunc (which
%   must be a function handle). Starts with initial solution 'initial' (a matrix). 
%   The values of each entry of the solution will be bounded by
%   the lower and upper matrices of the same dimensions.
%   Terminates when 5 successive solutions differ by
%   less than threshold.
%   temperature > 0 denotes the starting temperature of the function, which is
%   then mapped to coolingRate * temperature in successive iterations (ie.
%   nextTemp must be a strictly decreasing function for positive temperature.
%   numIt is the amount of times to attempt a pertubation per variable, ie.
%   for each temperature, and for variable, Rmax pertubations will be
%   attempted.
%   Based on approach described by Corana et al. (1987), Wong et al.
%   (2002), and Ogden et al. (2010)
%    coolingRate = 0.90;  % Recommended
    iterPerTemp = max(100, 5 * numel(initial));
    iterPerStep = 20;
    varyCrit = 2;
    % These could really be function parameters...
    % Using recommended values from Corana et al.
    
    assert(all(size(initial) == size(lower)) && ...
        all(size(lower) == size(upper)), ...
        'dimensions of initial, lower, and upper must be the same');
    for ii = 1:numel(initial)
        if lower(ii) > initial(ii)
            initial(ii) = lower(ii);
            warning('Initial value exceeds lower bound!! Forcing into bounds...');
        elseif upper(ii) < initial(ii)
            initial(ii) = upper(ii);
            warning('Initial value exceeds upper bound!! Forcing into bounds...');
        end
    end
    
    step = (upper - lower)/2;
    changeCount = zeros(size(initial));
    solHist = { initial initial initial initial };
    % Need a more effective way to force a first loop!!
    current = initial;
    phi = objFunc(current);
    optSoFar = current;
    optPhi = phi;
    iterCount = 0;
    while ((abs(phi - objFunc(solHist{1})) > threshold || ...
        abs(phi - objFunc(solHist{2})) > threshold) || ...
        (abs(phi - objFunc(solHist{3})) > threshold || ...
        abs(phi - objFunc(solHist{4})) > threshold)) || ...
        (iterCount < 5)
        % Tests if solutions differ by less than threshold.
        
        
        for stepInd = 1:iterPerTemp
            for iterInd = 1:iterPerStep
                candidate = current; 
                % Makes a copy.
                for varInd = 1:numel(initial)
                    candidate(varInd) = current(varInd) + ...
                         (2 * rand - 1) * (step(varInd));
                    while (candidate(varInd) < lower(varInd) || ...
                            candidate(varInd) > upper(varInd))
                        candidate(varInd) = current(varInd) + ...
                            (2 * rand - 1) * (step(varInd));
                    end
                    candidatePhi = objFunc(candidate);
                    deltaPhi = candidatePhi - objFunc(current);
                    if (deltaPhi <= 0 || exp(-deltaPhi / temperature) > rand) 
                        %evaluates second part as necessary
                        current = candidate;
                        phi = candidatePhi;
                        changeCount(varInd) = changeCount(varInd) + 1;
                        if (phi < optPhi)
                            % Keeping track of best solution so far...
                            optSoFar = current;
                            optPhi = phi;
                        end
                    end

                end
                
                
            end % end of iterPerStep
            % Update step size using formula from Corana et al
            for varInd = 1:numel(initial)
                if (changeCount(varInd) > 0.6 * iterPerStep)
                    step(varInd) = step(varInd) * ...
                        (1 + varyCrit * ...
                        ( changeCount(varInd) / iterPerStep - 0.6) / 0.4);
                elseif (changeCount(varInd) < 0.4 * iterPerStep)
                    step(varInd) = step(varInd) / ...
                        (1 + varyCrit * ...
                        (0.4 - changeCount(varInd) / iterPerStep) / 0.4);
                end
                
                %Check if step is absurdly large
                if (step(varInd) > upper(varInd) - lower(varInd))
                    step(varInd) = upper(varInd) - lower(varInd);
                end
            end
            changeCount = zeros(size(initial));
        end
        temperature = temperature * coolingRate;
        solHist = circshift(solHist, [0, 1]);
        solHist{1} = current;
        iterCount = iterCount + 1;
    end
    
    %display(solHist)
end

