clear variables;
format short g;
format compact;

ideal = 0;
regions = [3 4 7 8 9 10];

startingFrame = 11;
lastFrame = 28;
%Tac simulated from Pablo model
fullTAC = dlmread('Data/TACs/pabloModel_0vb/0noise/fullTAC0sigma.tac', '\t', 1, 0);
startTimes = fullTAC(:, 1);
endTimes = fullTAC(:, 2);
relTAC = fullTAC(:, regions); %the 6 brain regions
relTACint = cumtrapz(startTimes, relTAC);
relTAC = relTAC(startingFrame:end, :);
relTACint = relTACint(startingFrame:end, :);

sourceCp = dlmread('Data/Cps/pabloModel/pabloModel_0sigma.smpl', '\t', 1, 0);
bloodDrawFrame = 20; % 15th frame
startBlood = fullTAC(bloodDrawFrame, 1);
endBlood = fullTAC(bloodDrawFrame, 2);
bloodDrawTime = (startBlood  + endBlood) / 2;
singleBloodDraw = trapz(sourceCp(startBlood:endBlood, 2)) / ...
    (endBlood - startBlood);
bloodDrawErrFactor = 1e6;

ISAresult = ISA(relTAC, relTACint);

%Adjust for single blood draw
slope = (ISAresult(bloodDrawFrame - startingFrame + 1) - ...
    ISAresult(bloodDrawFrame - startingFrame)) / (endBlood - startBlood);
ISAresult = ISAresult * singleBloodDraw / slope;

iterResult = IterativeAlgorithm(relTAC, relTACint, ISAresult);

%Fit curve with biexponential model
times = (startTimes(startingFrame:end) + endTimes(startingFrame:end)) / 2;

[results, error] = fit2E(iterResult, times, singleBloodDraw, bloodDrawTime, ...
    bloodDrawErrFactor);

%%% Display data...
startTime = fullTAC(startingFrame, 1);
% sourceCp = sourceCp(sourceCp(:, 1) > startTime, :);
sourceCpInt = cumtrapz(sourceCp(:, 1), sourceCp(:, 2));
% Where to start intCp calculation?

CpIntAtTimes = interp1(sourceCp(:, 1), sourceCpInt, times);

Vt = calcVt(relTACint, relTAC, iterResult);

figure;
plot(times, CpIntAtTimes, ':', times, ISAresult, 'o-', ...
    times, iterResult, 'o-', times, model2E(results, times), 'o--');
legend('real Cp result', 'ISA result', 'It Alg result', 'biexponential model');

disp(['variables in model ' num2str(1) ':']);
disp(results(:, 1)');
disp(['Error in model' num2str(1) ' = ' num2str(error)]);
disp('Vt = ');
disp(Vt);

%combined graph
figure;
plot(times, CpIntAtTimes, ':',times, iterResult, 'bo', times, model2E(results, times), 'r-');
legend('real Cp int', 'generated Cp int', 'biexponential fit');

uisave;
