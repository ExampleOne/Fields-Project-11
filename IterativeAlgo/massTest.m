clear variables;
format short g;
format compact;

% sigmas = -2:0.2:2;
sigmas = -2;
% noises = 1:10;
noises = 8;
trials = 8;

numRegions = 7;
numNoises = length(noises);
% numTrials = 20;
numTrials = length(trials);

Vts = zeros(numRegions, length(sigmas), numNoises, numTrials);
DVR_table = zeros(numRegions - 1, length(sigmas), numNoises, numTrials);

for noiseInd = 1:numNoises
    for sigmaInd = 1:length(sigmas)
        tic;
        for trialInd = 1:numTrials
            
            sigma = sigmas(sigmaInd);
            TACPath = ['/Users/shixiaoying/Desktop/MATLAB/DVR/' ...
                'pabloModelTACs_cereb_ref_0vb/' num2str(noises(noiseInd))...
                'noise/' num2str(sigma) 'sigma/fullTACsf_' ...
                num2str(noises(noiseInd)) '_sim_' num2str(trials(trialInd)) '.tac'];
            CpPath = ['/Users/shixiaoying/Desktop/MATLAB/newtacs/Fields-Project-11/Data/Cps/pabloModel/pabloModel_' num2str(sigma) 'sigma.smpl'];
            Vts(:, sigmaInd, noiseInd, trialInd) = nihms(TACPath, CpPath, false, false);

            % Store the DVRs for this specific sigma and noise
            DVR_table(:, sigmaInd, noiseInd, trialInd) = ...
                Vts(1:numRegions - 1, sigmaInd, noiseInd, trialInd) ./  ...
                Vts(end, sigmaInd, noiseInd, trialInd);
        end
        toc;
    end
end

DVR_means = mean(DVR_table, 4);
DVR_stds = std(DVR_table, 0, 4);
DVR_maxes = max(DVR_table, [], 4);
DVR_mins = min(DVR_table, [], 4);


uisave;

