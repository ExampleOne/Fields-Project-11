clear variables;
format short g;
format compact;

startingFrame = 11;
lastFrame = 28;
%Tac simulated from Pablo model
fullTAC = dlmread('Data/TACs/pabloModel/0noise/fullTAC0sigma.tac', '\t', 1, 0);
trimmedTAC = fullTAC(startingFrame:end, :);
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

n = size(relTAC,2);
k = 1; %counter
diff = zeros(lastFrame - startingFrame + 1, 6, 6);
Cpint = zeros(lastFrame - startingFrame + 1, 6, 6);
dim = n*(n-1)/2;
Cpintlist = zeros(18,dim); %same thing with Cpint but in 2-dim
%va = zeros(4,dim);      %the [a1i,a2i,a1j,a2j] in Fields

%ISA part

for i = 1:n
    cti = relTAC(:,i);
   
    intcti = cumtrapz(startTimes,relTAC(:,i));
    
    %cpi = pdata(:,2);
    %cpi = cumtrapz(pdata(:,1),pdata(:,2));
    for j = i + 1:n
         ctj = relTAC(:,j);
         intctj = cumtrapz(startTimes,relTAC(:,j));
         % find the right singular vector va corresponding to smallest
         % right singular value of c4
         c4 = [cti, intcti, -ctj, -intctj];
         [U,S,V] = svd(c4);
         %va(:,k) = V(:,end);
         
         %evaluate the first guesses of Cpint
         diff(:,i,j) = abs(V(1,end)*cti + V(2,end)*intcti - (V(3,end)*ctj + V(4,end)*intctj));
         Cpint(:,i,j) = abs(V(1,end))*cti +abs(V(2,end))*intcti;
         Cpintlist(:,k) = Cpint(:,i,j);
         k = k+1;
    end
    
    

       
end

%find one cpint which minimize the summation of distances
dist = zeros(1,dim);

for l = 1:dim
    for m = 1:dim
        disttemp = norm(Cpintlist(:,l)-Cpintlist(:,m));
        dist(l) = dist(l) + disttemp;
    end
end

[M,minl] = min(dist);
Cpint1 = Cpintlist(:,minl);

%Adjust for single blood draw
slope = (Cpint1(bloodDrawFrame - startingFrame + 1) - ...
    Cpint1(bloodDrawFrame - startingFrame)) / (endTime - startTime);
Cpint1 = Cpint1 * singleBloodDraw / slope;

ISAresult = Cpint1(:);
va1 = zeros(2,n);
Cpint2 = zeros(18,n);
%%%%%begin the iterative algorithm part%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_err = 1;
err = 200;

while err > max_err
%regression to get new set of Vt and b
        for i = 1:n
            cti = relTAC(:,i);
            intcti = cumtrapz(startTimes,relTAC(:,i));
            Cpint1lm = fitglm([cti,intcti],Cpint1,'linear');
            vatemp = Cpint1lm.Coefficients.Estimate;
            va1(:,i) = vatemp(2:end);
            Cpint2(:,i) = va1(1,i)*cti + va1(2,i)*intcti;
        end   

        %take averages to estimate the cpint in new iteration----no
        %va2 = [mean(va1(1,:)),mean(va1(2,:))];

        %find the distance in loop
        dist = zeros(1,n);

        for l = 1:n
            for m = 1:n
                disttemp = norm(Cpint2(:,l)-Cpint2(:,m));
                dist(l) = dist(l) + disttemp;
            end
        end

        [M,minl] = min(dist);
        err = norm(Cpint1 - Cpint2(:,minl));
        Cpint1 = Cpint2(:,minl);
        va1 = zeros(2,n);
        Cpint2 = zeros(18,n);
end

% Extraplolate
% Cpint1ep = interp1(times, Cpint1, ox, 'linear', 'extrap');

%%use Simulated Annealing to fit the CPint1%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
times = (trimmedTAC(:,1) + trimmedTAC(:, 2)) / 2;

fity = @(b,x) b(1)*exp(-b(2)*x) + b(3)*exp(-b(4)*x)+b(5);


phi1 = @(b) trapz((Cpint1 - fity(b,times)).^2) + singleBloodDrawErrFactor * ...
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
CpIntAtTimes = CpIntAtTimes - CpIntAtTimes(1) + Cpint1(1);

figure;
plot(times, CpIntAtTimes, ':', times, ISAresult, 'o-', ...
    times, Cpint1, 'o-', times, fity(results, times), 'o--');
legend('real Cp result', 'ISA result', 'It Alg result', 'biexponential model');

error = phi1(results);%/trapz(Cpint1.^2) ;

% Includes integral of Y squared.

disp(['variables in model ' num2str(1) ':']);
disp(results(:, 1)');
disp(['Error in model' num2str(1) ' = ' num2str(error)]);
    
results2 = zeros(5,1); % change #variable
%real Cp data

%use model fitting ISA-italgo as fix point
shift = - CpIntAtTimes(1) + Cpint1(1);
%combined graph
figure;
plot(times, CpIntAtTimes, ':',times, Cpint1, 'bo', times, fity(results, times), 'r-');
legend('real Cp int', 'generated Cp int', 'fit generated Cp int');

uisave;
