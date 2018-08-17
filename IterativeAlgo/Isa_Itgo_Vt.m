clear variable;
%Tac simulated from Pablo model
fulltac = dlmread(['Data/TACs/pabloModel/0noise0vb/fullTAC0sigma.tac'], '\t', 1, 0);
startingFrame = 11;
lastFrame = 28;
plotpoints = lastFrame - startingFrame + 1;
totalpoints = size(fulltac,1);

trimmedTAC = fulltac(end - plotpoints + 1:end, :);
% pdata = plasmafile
datarelev = trimmedTAC(:,[3 4 7 8 9 10]); %the 6 brain regions
fulldatarelev =fulltac(:,[3 4 7 8 9 10]);

n = size(datarelev,2);
k = 1; %counter
diff = zeros(plotpoints,n,n);
Cpint = zeros(plotpoints,n,n);
dim = n*(n-1)/2;
Cpintlist = zeros(plotpoints,dim); %same thing with Cpint but in 2-dim
%va = zeros(4,dim);      %the [a1i,a2i,a1j,a2j] in Fields

intfulltac = zeros(totalpoints,n);

for i = 1:n
    intfulltac(:,i) = cumtrapz(fulltac(:,1),fulldatarelev(:,i));
    temp = intfulltac(:,i);
end

trimmedintTAC = intfulltac(end - plotpoints + 1:end,:);


sourceCp = dlmread('Data/Cps/pabloModel/pabloModel_0sigma.smpl', '\t', 1, 0);
bloodDrawFrame = 20; % 20th frame
startTime = fulltac(bloodDrawFrame, 1);
endTime = fulltac(bloodDrawFrame, 2);
bloodDrawTime = (startTime  + endTime) / 2;
singleBloodDraw = trapz(sourceCp(startTime:endTime, 2)) / ...
    (endTime - startTime);
singleBloodDrawErrFactor = 10000;

startTimes = trimmedTAC(:, 1);
endTimes = trimmedTAC(:, 2);


%ISA part
for i = 1:n
    cti = datarelev(:,i);
    intcti = trimmedintTAC(:,i);
    
    for j = i + 1 : n
         ctj = datarelev(:,j);
         intctj = trimmedintTAC(:,j);
         % find the right singular vector va corresponding to smallest
         % right singular value of c4
         c4 = [cti(2:end),intcti(2:end),ctj(2:end),intctj(2:end)];
         [U,S,V] = svd(c4);
         %va(:,k) = V(:,end);
         
         %evaluate the first guesses of Cpint
         diff(:,i,j) = abs(V(1,end)*cti + V(2,end)*intcti + (V(3,end)*ctj + V(4,end)*intctj));
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

%include the one sample blood data
slope = (Cpint1(bloodDrawFrame - startingFrame + 1) - ...
    Cpint1(bloodDrawFrame - startingFrame)) / (endTime - startTime);
Cpint1 = Cpint1 * singleBloodDraw / slope +725*60;
ISAresult = Cpint1(:)+725*60;


VtISA = zeros(n,1);
b = zeros(n,1);
for i = 1:n
        cti = datarelev(:,i);
        intcti = trimmedintTAC(:,i);
       
        %use logan plot expression to calculate Vt by lm
        dependentvariable = intcti(2:end)./cti(2:end);
        regressor = (Cpint1(2:end,:))./cti(2:end);
        [ERR,P] = fit_2D_data(regressor,dependentvariable,'no');
        VtISA(i,1) = P(1);
        b(i,1) = P(2);
end   

%PreditedTAC =

%%%%%begin the iterative algorithm part%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
va1 = zeros(2,n);
Cpint2 = zeros(18,n);
max_err = 1e-10;
err = 200;
iteration = 1;

%regression to get new set of Vt and b
while err > max_err

        
        for i = 1:n
            cti = datarelev(:,i);
            intcti = trimmedintTAC(:,i);
            Cpint1lm = fitglm([cti,intcti],Cpint1,'linear');
            vatemp = Cpint1lm.Coefficients.Estimate;
            va1(:,i) = vatemp(2:end);
            Cpint2(:,i) = va1(1,i)*cti + va1(2,i)*intcti;
           
        end   


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
        %va1 = zeros(2,n);
        %Cpint2 = zeros(18,n);
        iteration = iteration + 1;
end



%%use Simulated Annealing to fit the CPint1%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g;
times = trimmedTAC(:,1);

fity = @(b,x) b(1)*(1-exp(-b(2)*x)) + b(3)*(1-exp(-b(4)*x))+b(5);

errors = zeros(1,2);
results = zeros(5,2); %the 4 can be changed 

phi1 = @(b) trapz((Cpint1 - fity(b,times)).^2);

tic;
results(:,1) = simulatedAnnealing(phi1,...
    [0.1 0.00001 0.1 0.00001 0], ...
    [-1e9 0 -1e9 0 -1e9], ...
    [1e9 0.05 1e9 0.05 1e9 ], 1000, 1e-7, 0.85);
disp(['Time elapsed in round ' num2str(1) ':']);
toc;

datap = dlmread(['Data/Cps/pabloModel/pabloModel_0sigma.smpl'], '\t', 1, 0);
otp = datap(:,1);
oyp = datap(:,2);
oypint = cumtrapz(otp,oyp);

yp = oyp(otp >= 1500);
tp = otp(otp >= 1500);
ypint = cumtrapz(tp, yp);

CpIntAtTimes = interp1(oypint, times);

%real Cp data
%%%%%%%%%%%%%%%%%%%calculate the factor for Cpint%%%%%%%%%%%%%
fitCpintAtSample = fity(results(:, 1), endTime)-fity(results(:, 1), startTime);
Cpintfactor = trapz(sourceCp(startTime:endTime, 2))/fitCpintAtSample;

%%%%%%%%%%%%%%%%%%%calculate the factor for Cp%%%%%%%%%%%%%%%%
fitCpAtSample = results(1,1)*results(2,1)*exp(-results(2,1)*bloodDrawTime)+...
                results(3,1)*results(4,1)*exp(-results(4,1)*bloodDrawTime);
Cpfactor = singleBloodDraw/fitCpAtSample;



% fitp = @(b,x) b(1)*exp(-b(2)*x) + b(3)*exp(-b(4)*x)+b(5); %#change format
% 
% phi2 = @(b) trapz((ypint - fitp(b,tp)).^2) ;%- trapz(;
% 
% tic;
% results2(:,1) = simulatedAnnealing(phi2,...
%    [0.1 0.00001 0.1 0.00001 0], ...
%     [-50000 0 -50000 0 -50000 ], ...
%     [1000 0.001 1000 0.001 50000 ], 1000, 1e-7, 0.85);
% disp(['Time elapsed in round ' num2str(2) ':']);
% toc;
% 
% figure;
% plot(otp, oyp, ':', tp, fitp(results2(:, 1), tp), '--');
% legend(['data 2' ], ['model 2' ]);
% 
% errors (1,2) = phi2(results2(:, 1))/trapz(ypint.^2) ;
% 
% 
% disp(['variables in model ' num2str(2) ':']);
% disp(results2(:, 1)');
% disp(['Error in model' num2str(2) ' = ' num2str(errors(1,2))]);

%use model fitting ISA-italgo as fix point
figure;
plot(otp, oypint, ':', times, ISAresult, 'o-', ...
    times, Cpint1, 'o-');
legend('real Cp result', 'ISA result', 'It Alg result 1');

%errors (1,1) = phi1(results(:, 1))/trapz(Cpint1.^2) ;

%Includes integral of Y squared.

disp(['variables in model ' num2str(1) ':']);
disp(results(:, 1)');
disp(['Error in model' num2str(1) ' = ' num2str(errors(1,1))]);
    


%shift3 = oypint(end)-Cpint1(end);
%combined graph
figure;
plot(otp, oypint, ':',times, Cpint1, 'bo', times, fity(results(:, 1), times), 'r-');
legend(['real Cp int' ], ['gened Cp int'], ['fit gened Cp int']);
% 

%%%%calculate the Vt 
Vt = zeros(n,1);
for i = 1:n
        cti = datarelev(:,i);
        intcti = cumtrapz(trimmedTAC(:,1),datarelev(:,i));
       
        %use logan plot expression to calculate Vt by lm
        dependentvariable = intcti(2:end)./cti(2:end);
        regressor = (Cpint1(2:end,:)*Cpintfactor)./cti(2:end);
        [ERR,P] = fit_2D_data(regressor,dependentvariable,'no');
        coeffest = P(1);
        Vt(i,1) = coeffest;
end   

disp('Vt:=');
disp(Vt);
disp('Vt from ISA = ');
disp(VtISA);

uisave;