clear variables;
format short g;
format compact;

sigmas = -2:0.2:2;
noises = 0:10;

googleSpreadSheetID = '1mBe8fQp-5RCBkfMJ35gHHYcZTxHu5e9TeM7xo15F85o';
googleSheetID = '1200493859';

numRegions = 6;
numNoises = length(noises);
numTrials = 200;
numSigmas = length(sigmas);

Vts = zeros(numRegions, length(sigmas), numNoises, numTrials);
DVRs = zeros(numRegions - 1, length(sigmas), numNoises, numTrials);

for noiseInd = 1:numNoises
    for sigmaInd = 1:numSigmas
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

meansOut = reshape(VtMeans(:, :, :), [numRegions * numSigmas numNoises]);
stdsOut = reshape(VtStds(:, :, :), [numRegions * numSigmas numNoises]);
maxesOut = reshape(VtMaxes(:, :, :), [numRegions * numSigmas numNoises]);
minsOut = reshape(VtMins(:, :, :), [numRegions * numSigmas numNoises]);

mat2sheets(googleSpreadSheetID, googleSheetID, [13 4], meansOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 18], stdsOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 32], minsOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 46], maxesOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [1 3], {['Last update: ' datestr(now)]});

uisave;

