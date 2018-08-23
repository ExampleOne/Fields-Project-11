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
time = zeros(numNoises, length(sigmas));

for noiseInd = 1:numNoises
    for sigmaInd = 1:length(sigmas)
        tic;
        for trialInd = 1:numTrials
            sigma = sigmas(sigmaInd);
            TACPath = ['/home/jouyang/Downloads/GIT_code/'...
                'pabloModelTACs_cereb_ref_0vb_200_highK1/' num2str(noises(noiseInd))...
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
        time(noiseInd,sigmaInd) = toc;
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

googleSpreadSheetID = '147SxPtM9_ycQNbDtGV7IZ7BykBzEYpm3G8kPx2KdtUI';
googleSheetID = '1071545489';

meansOut = reshape(DVR_means(:, :, :), [(numRegions-1) * length(sigmas) numNoises]);
stdsOut = reshape(DVR_stds(:, :, :), [(numRegions-1) * length(sigmas) numNoises]);
minsOut = reshape(DVR_mins(:, :, :), [(numRegions-1) * length(sigmas) numNoises]);
maxesOut = reshape(DVR_maxes(:, :, :), [(numRegions-1) * length(sigmas) numNoises]);

mat2sheets(googleSpreadSheetID, googleSheetID, [13 4], meansOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 18], stdsOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 32], minsOut);
mat2sheets(googleSpreadSheetID, googleSheetID, [13 46], maxesOut);



Totaltime = sum(sum(time));
disp(['Total time = ' num2str(Totaltime)]);

for noiseInd = 1:numNoises
    disp(['noiseInd = ' num2str(noiseInd)]);
    for sigmaInd = 1:length(sigmas)
        if time(noiseInd,sigmaInd) > 20
            disp(['noise index = ' num2str(noiseInd)]);
            disp(['sigma index = ' num2str(sigmaInd)]);
            disp(['abnormal time = ' num2str(time(noiseInd,sigmaInd))]);
        end
    end
end

uisave;

