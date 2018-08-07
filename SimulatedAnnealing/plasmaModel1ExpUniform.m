% Fitting an exponential model to the plasma data.

data = dlmread('lammerplasma.smpl', '\t', 1, 0);

trimmed_data = trimmed(68:end,:);



% Single exponential model
singleExp = @(c, t) c(1) * exp(-c(2) * t) + c(3);
phiSE = @(c) norm(singleExp(c, trimmed_data(:, 1)) - trimmed_data(:, 2));
tic;
resultsSE = simulatedAnnealing(phiSE, [5 0.01 5], ...
    [-1e4 0 0], ...
    [1e4 0.05 1000], 1e7, 1e-5, 0.9);
toc;
display(resultsSE);
display(phiSE(resultsSE));

figure;
plot(data(:,1), data(:, 2), ':', data(:,1), singleExp(resultsSE, data(:, 1)), '--');
legend('data', 'model');





