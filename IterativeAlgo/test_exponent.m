clear variables;
format short g;
format compact;

startingFrame = 11;
lastFrame = 28;
%Tac simulated from Pablo model
fullTAC = dlmread('Data/TACs/pabloModel/0noise/fullTAC0sigma.tac', '\t', 1, 0);
trimmedTAC = fullTAC(startingFrame:lastFrame, :);
% pdata = plasmafile
relTAC = trimmedTAC(:,[3 4 7 8 9 10]); %the 6 brain regions

sourceCp = dlmread('Data/Cps/pabloModel/pabloModel_0sigma.smpl', '\t', 1, 0);
bloodDrawFrame = 20; % 15th frame
startTime = fullTAC(bloodDrawFrame, 1);
endTime = fullTAC(bloodDrawFrame, 2);
bloodDrawTime = (startTime  + endTime) / 2;
singleBloodDraw = trapz(sourceCp(startTime:endTime, 2)) / ...
    (endTime - startTime);
bloodDrawErrFactor = 1e6;

startTimes = trimmedTAC(:, 1);
endTimes = trimmedTAC(:, 2);

ISAresult = ISA(relTAC, startTimes);

%Adjust for single blood draw
slope = (ISAresult(bloodDrawFrame - startingFrame + 1) - ...
    ISAresult(bloodDrawFrame - startingFrame)) / (endTime - startTime);
ISAresult = ISAresult * singleBloodDraw / slope;

iterResult = IterativeAlgorithm(relTAC, ISAresult, startTimes);

%%use Simulated Annealing to fit the CPint1%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
times = (trimmedTAC(:,1) + trimmedTAC(:, 2)) / 2;

results = fit2E(iterResult, times, singleBloodDraw, bloodDrawTime, ...
    bloodDrawErrFactor);

%%% Display data...
sourceCpInt = cumtrapz(sourceCp(:, 1), sourceCp(:, 2));

CpIntAtTimes = interp1(sourceCpInt, times);
CpIntAtTimes = CpIntAtTimes - CpIntAtTimes(1) + iterResult(1);

figure;
plot(times, CpIntAtTimes, ':', times, ISAresult, 'o-', ...
    times, iterResult, 'o-', times, model2E(results, times), 'o--');
legend('real Cp result', 'ISA result', 'It Alg result', 'biexponential model');

error = phi1(results);%/trapz(Cpint1.^2) ;

% Includes integral of Y squared.

disp(['variables in model ' num2str(1) ':']);
disp(results(:, 1)');
disp(['Error in model' num2str(1) ' = ' num2str(error)]);

%combined graph
figure;
plot(times, CpIntAtTimes, ':',times, iterResult, 'bo', times, model2E(results, times), 'r-');
legend('real Cp int', 'generated Cp int', 'biexponential fit');

uisave;
