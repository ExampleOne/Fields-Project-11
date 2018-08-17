function [CpInt, Vt] = IterativeAlgorithm(TAC, TACint, initialCp)
%ITERATIVEALGORITHM Performs the Iterative Algorithm described in NIHMS
%   Repetetively performs a linear regression to predict CpInt
    numRegions = size(TAC, 2);
    Vts = zeros(numRegions, 5);
    bs = zeros(numRegions, 5);
    threshold = 1e-4;
    
    shouldContinue = true;
    CpInt = initialCp;
    CpEstimates = zeros(length(CpInt), numRegions);
    deltas = zeros(1, numRegions);
    errs = zeros(1, numRegions);
    while shouldContinue
        %Calculate Vt and b
        for ii = 1:numRegions
            [err, B] = fit_2D_data(CpInt ./ TAC(:, ii), ...
                TACint(:, ii) ./ TAC(:, ii), 'no');
            Vts(ii, 1) = B(1);
            bs(ii, 1) = B(2);
            errs(ii) = err;
            CpEstimates(:, ii) = -bs(ii, 1) / Vts(ii, 1) * TAC(:, ii) + ...
                1 / Vts(ii, 1) * TACint(:, ii);
        end
        for ii = 1:numRegions
            deltas(ii) = norm(repmat(CpEstimates(:, ii), 1, numRegions) ...
                - CpEstimates, 'fro'); 
            % Select best Cp on the one closest to the others...
            % Consider other approaches?
        end
        [~, minIndex] = min(deltas);
        CpInt = CpEstimates(:, minIndex);
        
        shouldContinue = norm(Vts - circshift(Vts, [0 1]), 'fro')^2 + ...
            norm(bs - circshift(bs, [0 1]), 'fro')^2 >= threshold;
        Vts = circshift(Vts, [0 1]);
        bs = circshift(bs, [0 1]);
    end
    
    Vt = Vts(:, 1);
    % CpInt already at its best.
end

