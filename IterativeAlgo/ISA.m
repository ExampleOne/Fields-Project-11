function bestCpInt = ISA(TACs, startTimes)
%ISA Performs the Intersectional Search Algorithm described in NIHMS
%   Finds the intersection of two planes to estimate Cp.
    n = size(TACs, 2);
    k = 1; %counter
    Cpint = zeros(size(TACs, 1), size(TACs, 2), size(TACs, 2));
    dim = n*(n-1)/2;
    Cpintlist = zeros(18,dim); %same thing with Cpint but in 2-dim
    %va = zeros(4,dim);      %the [a1i,a2i,a1j,a2j] in Fields

    for i = 1:n
        cti = TACs(:,i);

        intcti = cumtrapz(startTimes, TACs(:,i));

        %cpi = pdata(:,2);
        %cpi = cumtrapz(pdata(:,1),pdata(:,2));
        for j = i + 1:n
             ctj = TACs(:,j);
             intctj = cumtrapz(startTimes,TACs(:,j));
             % find the right singular vector va corresponding to smallest
             % right singular value of c4
             c4 = [cti, intcti, -ctj, -intctj];
             [~, ~, V] = svd(c4);
             %va(:,k) = V(:,end);

             %evaluate the first guesses of Cpint
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

    [~, minl] = min(dist);
    bestCpInt = Cpintlist(:,minl);
end

