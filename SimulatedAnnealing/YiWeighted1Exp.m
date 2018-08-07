% Fitting an exponential model to the plasma data.

clear variables;

data = dlmread('ModellingCp/plasmaFunctions/plasma1.smpl', '\t', 1, 0);

% Single exponential model
singleExp = @(c, t) c(1) * exp(-c(2) * t) + c(3);
phi = @(c) YiWeightedPhi(data, @(t) singleExp(c, t));

tic;
results = simulatedAnnealing(...
    phi, ...
    [5 0.01 5], ...
    [-1e4 0 0], ...
    [1e4 0.05 1000], 1e7, 1e-5, 0.9);
toc;
display(results);
display(phi(results));

figure;
plot(data(:,1), data(:, 2), data(:,1), singleExp(results, data(:, 1)));
legend('data', 'model');





