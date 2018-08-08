function Vt = calcVt(TACint, TACs, Cpint)
%CALCVT Calculates Vt
%   Uses Logan method.
    Vt = zeros(size(TACs, 2), 1);
    for ii = 1:size(TACs, 2)
        X = [ones(size(TACs, 1), 1) Cpint ./ TACs(:, ii)];
        B = X \ (TACint(:, ii) ./ TACs(:, ii));
        Vt(ii) = B(2);
    end
end

