clear variables;

data = dlmread('lammerplasma.smpl', '\t', 1, 0);

trimmed_data = trimmed(68:end,:);
%Double exponential model

doubleExp = @(c, t) c(1) * exp(-c(2) * t) + c(3) * exp(-c(4) * t) + ...
    c(5);
phiDE = @(c) norm(doubleExp(c, trimmed_data(:, 1)) - trimmed_data(:, 2));
tic;
resultsDE = simulatedAnnealing(phiDE, [5 0.01 5 0.01 5], ...
    [-1e4 0 -1e4 0 0], [1e4 0.05 1e4 0.05 1000], 1e7, 1e-5, 0.9);
toc;
display(resultsDE);
display(phiDE(resultsDE));

figure
figure;
plot(data(:,1), data(:, 2), ':', data(:,1), doubleExp(resultsDE, data(:, 1)), '--');
legend('data', 'model');