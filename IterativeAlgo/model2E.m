function Cp = model2E(b, t)
%MODEL2E The biexponential model used to fit the curve.
%   Detailed explanation goes here
    Cp = b(1) * exp(-b(2) * t) + b(3) * exp(-b(4)*t) + b(5);
end

