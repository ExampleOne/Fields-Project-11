clear variables;
format short g;
format compact;

sigmas = -2:0.2:2;
noises = 0:10;

numRegions = 6;
numNoises = length(noises);
numTrials = 200;

Vts = zeros(numRegions, length(sigmas), numNoises, numTrials);
DVRs = zeros(numRegions - 1, length(sigmas), numNoises, numTrials);

for noiseInd = 1:numNoises
    disp(['Looking at noise level: ' num2str(noises(noiseInd))]);
    for sigmaInd = 1:length(sigmas)
        tic;
        for trialInd = 1:numTrials
            sigma = sigmas(sigmaInd);
            TACPath = ['/home/qtupker/Documents/Fields Project 11/' ...
                'workspace4/Square3EinputFunction/' ...
                'pabloModelTACs_cereb_ref_0vb_with_negatives/' num2str(noises(noiseInd))...
                'noise/' num2str(sigma) 'sigma/fullTACsf_' ...
                num2str(noises(noiseInd)) '_sim_' num2str(trialInd) '.tac'];
            CpPath = ['Data/Cps/pabloModel/pabloModel_' num2str(sigma) 'sigma.smpl'];
            Vts(:, sigmaInd, noiseInd, trialInd) = nihms(TACPath, CpPath, false, false);
        end
        toc;
    end
end

VtMeans = mean(Vts, 4);
VtStds = std(Vts, 0, 4);
VtMaxes = max(Vts, [], 4);
VtMins = min(Vts, [], 4);

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(VtMeans(:, sigmaInd, noiseInd));
    end
end

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(VtStds(:, sigmaInd, noiseInd));
    end
end

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(VtMins(:, sigmaInd, noiseInd));
    end
end

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        disp(VtMaxes(:, sigmaInd, noiseInd));
    end
end

uisave;

