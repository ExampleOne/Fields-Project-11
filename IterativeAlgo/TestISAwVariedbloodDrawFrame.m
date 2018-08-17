clear variable;
%Tac simulated from Pablo model
fulltac = dlmread(['Data/TACs/pabloModel/0noise0vb/fullTAC-2sigma.tac'], '\t', 1, 0);
sourceCp = dlmread('Data/Cps/pabloModel/pabloModel_0.2sigma.smpl', '\t', 1, 0);

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


frame = startingFrame:lastFrame;
lframe = size(frame,2);
VtISA = zeros(6,lframe-2);

for i = startingFrame + 1:1:lastFrame-1
    
    bloodDrawFrame = i; % 20th frame
    
    startTime = fulltac(bloodDrawFrame, 1);
    endTime = fulltac(bloodDrawFrame, 2);
    bloodDrawTime = (startTime  + endTime) / 2;
    singleBloodDraw = trapz(sourceCp(startTime:endTime, 2)) / ...
        (endTime - startTime);
    singleBloodDrawErrFactor = 10000;

    startTimes = trimmedTAC(:, 1);
    endTimes = trimmedTAC(:, 2);


    %ISA part
    for m = 1:n
        cti = datarelev(:,m);
        intcti = trimmedintTAC(:,m);

        for j = m + 1 : n
             ctj = datarelev(:,j);
             intctj = trimmedintTAC(:,j);
             % find the right singular vector va corresponding to smallest
             % right singular value of c4
             c4 = [cti(2:end),intcti(2:end),ctj(2:end),intctj(2:end)];
             [U,S,V] = svd(c4);
             %va(:,k) = V(:,end);

             %evaluate the first guesses of Cpint
             diff(:,m,j) = abs(V(1,end)*cti + V(2,end)*intcti + (V(3,end)*ctj + V(4,end)*intctj));
             Cpint(:,m,j) = abs(V(1,end))*cti +abs(V(2,end))*intcti;
             Cpintlist(:,k) = Cpint(:,m,j);
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

    for k = 1:n
            cti = datarelev(:,k);
            intcti = trimmedintTAC(:,k);

            %use logan plot expression to calculate Vt by lm
            dependentvariable = intcti(2:end)./cti(2:end);
            regressor = (Cpint1(2:end,:))./cti(2:end);
            [ERR,P] = fit_2D_data(regressor,dependentvariable,'no');
            VtISA(k,i) = P(1);
            b(k,i) = P(2);
    end   

end

 disp('Vt from ISA = ');
 disp(VtISA);