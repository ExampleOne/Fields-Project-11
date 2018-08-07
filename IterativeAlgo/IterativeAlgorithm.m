function bestCpInt = IterativeAlgorithm(TACs, initialCp, startTimes)
%ITERATIVEALGORITHM Performs the Iterative Algorithm described in NIHMS
%   Repetetively performs a linear regression to predict CpInt
    numRegions = size(TACs, 2);
    va1 = zeros(2, numRegions);
    Cpint2 = zeros(18, numRegions);
    max_err = 1;
    err = 200;
    bestCpInt = initialCp;

    while err > max_err
    %regression to get new set of Vt and b
            for i = 1:numRegions
                cti = TACs(:,i);
                intcti = cumtrapz(startTimes, TACs(:,i));
                Cpint1lm = fitglm([cti,intcti], bestCpInt,'linear');
                vatemp = Cpint1lm.Coefficients.Estimate;
                va1(:,i) = vatemp(2:end);
                Cpint2(:,i) = va1(1,i)*cti + va1(2,i)*intcti;
            end   

            %take averages to estimate the cpint in new iteration----no
            %va2 = [mean(va1(1,:)),mean(va1(2,:))];

            %find the distance in loop
            dist = zeros(1,numRegions);

            for l = 1:numRegions
                for m = 1:numRegions
                    disttemp = norm(Cpint2(:,l)-Cpint2(:,m));
                    dist(l) = dist(l) + disttemp;
                end
            end

            [~, minl] = min(dist);
            err = norm(bestCpInt - Cpint2(:,minl));
            bestCpInt = Cpint2(:,minl);
            va1 = zeros(2,numRegions);
            Cpint2 = zeros(18,numRegions);
    end
end

