clear variables;
format short g;
format compact;

sigmas = -2:0.2:2;
noises = 0:10;

numRegions = 7;
numNoises = length(noises);
numTrials = 20;

Vts = zeros(numRegions, length(sigmas), numNoises, numTrials);
DVRs = zeros(numRegions - 1, length(sigmas), numNoises, numTrials);

for noiseInd = 1:numNoises
    for sigmaInd = 1:length(sigmas)
        tic;
        for trialInd = 1:numTrials
            sigma = sigmas(sigmaInd);
            TACPath = ['/home/qtupker/Documents/Fields Project 11/workspace4/Square3EinputFunction/' ...
                'pabloModelTACs_cereb_ref_0vb/' num2str(noises(noiseInd))...
                'noise/' num2str(sigma) 'sigma/fullTACsf_' ...
                num2str(noises(noiseInd)) '_sim_' num2str(trialInd) '.tac'];
            CpPath = ['Data/Cps/pabloModel/pabloModel_' num2str(sigma) 'sigma.smpl'];
            Vts(:, sigmaInd, noiseInd, trialInd) = nihms(TACPath, CpPath, false, false);

            % Store the DVRs for this specific sigma and noise
            DVRs(:, sigmaInd, noiseInd, trialInd) = ...
                Vts(1:numRegions - 1, sigmaInd, noiseInd, trialInd) ./  ...
                Vts(end, sigmaInd, noiseInd, trialInd);
        end
        toc;
    end
end

DVR_means = mean(DVRs, 4);
DVR_stds = std(DVRs, 0, 4);

uisave;

