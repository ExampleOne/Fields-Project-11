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
singleBloodDrawErrFactor = 10000;

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

fity = @(b,x) b(1)*exp(-b(2)*x) + b(3)*exp(-b(4)*x)+b(5);


phi1 = @(b) trapz((iterResult - fity(b,times)).^2) + singleBloodDrawErrFactor * ...
    ( (fity(b, bloodDrawTime + 0.05) - fity(b, bloodDrawTime - .05)) * 10 ...
    - singleBloodDraw) ^ 2;
% Takes single blood draw into account for error.

tic;
results = simulatedAnnealing(phi1,...
    [0.1 0.00001 0.1 0.00001 0], ...
    [-1e6 0 -1e6 0 0], ...
    [0 0.05 0 0.05 2e6 ], 1000, 1e-7, 0.85);
disp(['Time elapsed in round ' num2str(1) ':']);
toc;

otp = sourceCp(:,1);
oyp = sourceCp(:,2);
oypint = cumtrapz(otp,oyp);

yp = oyp(otp >= 1500);
tp = otp(otp >= 1500);
ypint = cumtrapz(tp, yp);

CpIntAtTimes = interp1(oypint, times);
CpIntAtTimes = CpIntAtTimes - CpIntAtTimes(1) + iterResult(1);

figure;
plot(times, CpIntAtTimes, ':', times, iterResult, 'o-', ...
    times, iterResult, 'o-', times, fity(results, times), 'o--');
legend('real Cp result', 'ISA result', 'It Alg result', 'biexponential model');

error = phi1(results);%/trapz(Cpint1.^2) ;

% Includes integral of Y squared.

disp(['variables in model ' num2str(1) ':']);
disp(results(:, 1)');
disp(['Error in model' num2str(1) ' = ' num2str(error)]);
    
results2 = zeros(5,1); % change #variable
%real Cp data

%use model fitting ISA-italgo as fix point
shift = - CpIntAtTimes(1) + iterResult(1);
%combined graph
figure;
plot(times, CpIntAtTimes, ':',times, iterResult, 'bo', times, fity(results, times), 'r-');
legend('real Cp int', 'generated Cp int', 'fit generated Cp int');

uisave;
