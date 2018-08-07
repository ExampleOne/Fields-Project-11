data = dlmread('ModellingCp/plasma1.smpl', '\t', 1, 0);

trimmed_data = data;%(68:end,:); Manually trimmed file...

% Quadruple exponential model
model4Exp = @(c, t) c(1) * exp(-c(2) * t) + ...
    c(3) * exp(-c(4) * t) + ...
    c(5) * exp(-c(6) * t) + ...
    c(7) * exp(-c(8)) + c(9);
phi4Exp = @(c) norm(model4Exp(c, trimmed_data(:, 1)) - trimmed_data(:, 2)); %plasmaPhi(data, model, c);
tic;
results = simulatedAnnealing(phi4Exp, ...
    [100 0.01 100 0.01 100 0.01 100 0.01 100], ...
    [-1e5 0 -1e5 0 -1e5 0 -1e5 0 0], ...
    [1e5 0.05 1e5 0.05 1e5 0.05 1e5 0.05 1000], 1e10, 1e-5, 0.95);
toc;
display(results);
display(phi4Exp(results));
figure;
plot(data(:,1), data(:, 2), ':', data(:,1), model4Exp(results, data(:, 1)), '--');
legend('data', 'model');