clear variables;
format short g;
format compact;

numRegions = 6;
numNoises = 11;
numTrials = 20;

sigmas = -2:0.2:2;
Vts = zeros(numRegions, numNoises, length(sigmas), numTrials);

tic;
for sigmaInd = 1:length(sigmas)
    
    sigma = sigmas(sigmaInd);
    TACPath = ['/home/qtupker/Documents/Fields Project 11/workspace4/' ...
        'Square3EinputFunction/pabloModelTACs_0vb/0noise/fullTAC' ...
        num2str(sigma) 'sigma.tac'];
    CpPath = ['/home/qtupker/Documents/Fields Project 11/workspace4/' ...
        'Square3EinputFunction/Cps/pabloModel_' num2str(sigma) 'sigma.smpl'];
    Vts(:, 1, sigmaInd, 1) = nihms(TACPath, CpPath, false, false);
    
end
toc;
    
    




uisave;

