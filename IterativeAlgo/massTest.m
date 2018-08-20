clear variables;
format short g;
format compact;

numRegions = 6;
numNoises = 11;
numTrials = 20;

sigmas = -2:0.2:2;
Vts = zeros(numRegions, numNoises, length(sigmas), numTrials);


for sigmaInd = 1:length(sigmas)
    tic;
    sigma = sigmas(sigmaInd);
    TACPath = ['/home/qtupker/Documents/Fields Project 11/workspace4/' ...
        'Square3EinputFunction/pabloModelTACs_0vb/0noise/fullTAC' ...
        num2str(sigma) 'sigma.tac'];
    CpPath = ['/home/qtupker/Documents/Fields Project 11/workspace4/' ...
        'Square3EinputFunction/Cps/pabloModel_' num2str(sigma) 'sigma.smpl'];
    Vts(:, 1, sigmaInd, 1) = nihms(TACPath, CpPath, true, false);
    toc;
end

Vts = permute(Vts, [1 3 2 4]);
Vts = Vts(:, :, 1, 1);
    
    




uisave;

