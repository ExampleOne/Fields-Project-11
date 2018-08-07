function phi = YiWeightedPhi(data, model)
%YIWIEGHTEDPHI Calculates the sum of square differences weighted by Yi.
%   That is, if Yi represents the plasma data, we calculate the difference
%   expressed by 
%   $$ \sum_i Y_i (Y_i - f(t_i))^2. $$
%   Here data is Nx2 matrix, with one column of times, and one column of
%   values (Yi), and model is the function f.
    phi = sum(data(:, 2) .* (data(:, 2) - model(data(:, 1))) .^ 2);
end

