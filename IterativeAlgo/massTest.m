clear variables;
format short g;
format compact;

sigmas = -2:0.2:2;
noises = 0:10;

numRegions = 7;
numNoises = length(noises);
numTrials = 200;

Vts = zeros(numRegions, length(sigmas), numNoises, numTrials);
DVR_table = zeros(numRegions - 1, length(sigmas), numNoises, numTrials);

for noiseInd = 1:numNoises
    for sigmaInd = 1:length(sigmas)
        tic;
        for trialInd = 1:numTrials
            sigma = sigmas(sigmaInd);
            TACPath = ['/home/jouyang/Downloads/GIT_code/'...
                'pabloModelTACs_cereb_ref_0vb_with_negatives/' num2str(noises(noiseInd))...
                'noise/' num2str(sigma) 'sigma/fullTACsf_' ...
                num2str(noises(noiseInd)) '_sim_' num2str(trialInd) '.tac'];
            CpPath = ['/home/jouyang/Downloads/GIT_code/'...
                'Data/Cps/pabloModel/pabloModel_' num2str(sigma) 'sigma.smpl'];
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

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(DVR_means(:, sigmaInd, noiseInd));
    end
end
% 
for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(DVR_stds(:, sigmaInd, noiseInd));
    end
end
% 
for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(DVR_mins(:, sigmaInd, noiseInd));
    end
end
% 
for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(DVR_maxes(:, sigmaInd, noiseInd));
    end
end



uisave;

