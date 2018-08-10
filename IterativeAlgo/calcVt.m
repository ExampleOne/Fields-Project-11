function Vt = calcVt(TACint, TACs, Cpint)
%CALCVT Calculates Vt
%   Uses Logan method.
    Vt = zeros(size(TACs, 2), 1);
    for ii = 1:size(TACs, 2)
        [~, B] = fit_2D_data(Cpint ./ TACs(:, ii), ...
            TACint(:, ii) ./ TACs(:, ii), 'no');
        Vt(ii) = B(1);
    end
end

