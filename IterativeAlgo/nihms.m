function Vt = nihms(TACpath, CpPath, showGraphs, showText)
%NIHMS implements the method described in the NIHMS paper
%

regions = [3 4 7 8 9 10];

startingFrame = 11;
fullTAC = dlmread(TACpath, '\t', 1, 0);
startTimes = fullTAC(:, 1);
endTimes = fullTAC(:, 2);
TAC = fullTAC(:, regions); %the 6 brain regions
TACint = cumsum((endTimes - startTimes) .* TAC); %integrate
TAC = TAC(startingFrame:end, :);
TACint = TACint(startingFrame:end, :);

sourceCp = dlmread(CpPath, '\t', 1, 0);
bloodDrawFrame = 20; % 20th frame
startBlood = startTimes(bloodDrawFrame);
endBlood = endTimes(bloodDrawFrame);
bloodDrawTime = (startBlood  + endBlood) / 2;
singleBloodDraw = trapz(sourceCp(startBlood:endBlood, 2)) / ...
    (endBlood - startBlood);
bloodDrawErrFactor = 1e6;

%CALCULATING AUC!!!
earlyIndices = sourceCp(:, 1) < startTimes(startingFrame);
auc = trapz(sourceCp(earlyIndices, 1), sourceCp(earlyIndices, 2));

ISAresult = ISA(TAC, TACint);

%Adjust for single blood draw
shortRegressRadius = 2;
startIndex = bloodDrawFrame - shortRegressRadius;
endIndex = bloodDrawFrame + shortRegressRadius + 1; % since we draw blood through a whole frame...
[~, P] = fit_2D_data(startTimes(startIndex:endIndex), ...
    ISAresult(startIndex - startingFrame + 1:endIndex - startingFrame + 1), 'no');
ISAslope = P(1);

ISAresult = ISAresult * singleBloodDraw / ISAslope;


ISAresult = ISAresult + auc - ISAresult(1);
% ie. we *fix* the first point of ISA!!! This is quite important.

[iterResult, Vt] = IterativeAlgorithm(TAC, TACint, ISAresult);

if showGraphs
    %Fit curve with biexponential model'
    CpTimes = (startTimes(6:end) + endTimes(6:end)) / 2;
    times = (startTimes(startingFrame:end) + endTimes(startingFrame:end)) / 2;

%     [results, error] = fit2E(iterResult, times, singleBloodDraw, bloodDrawTime, ...
%         bloodDrawErrFactor);

    
    %%% Display data...
    % sourceCp = sourceCp(sourceCp(:, 1) > startTime, :);
    sourceCpInt = cumtrapz(sourceCp(:, 1), sourceCp(:, 2));
    % Where to start intCp calculation?

    CpIntAtTimes = interp1(sourceCp(:, 1), sourceCpInt, CpTimes);
    
    figure;
    plot(CpTimes, CpIntAtTimes, ':', times, ISAresult, 'o-', ...
        times, iterResult, 'o-')%, times, model2E(results, times), 'o--');
    legend('real Cp result', 'ISA result', 'It Alg result', 'biexponential model');
    
    %combined graph
%     figure;
%     plot(times, CpIntAtTimes, ':',times, iterResult, 'bo')%, times, model2E(results, times), 'r-');
%     legend('real Cp int', 'generated Cp int', 'biexponential fit');
end

if showText
    disp(['variables in model ' num2str(1) ':']);
    disp(results(:, 1)');
    disp(['Error in model' num2str(1) ' = ' num2str(error)]);
    disp('Vt = ');
    disp(Vt);
end

end
