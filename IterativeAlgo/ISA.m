function bestCpInt = ISA(TAC, TACint)
%ISA Performs the Intersectional Search Algorithm described in NIHMS
%   Finds the intersection of two planes to estimate Cp.
    CpInts = [];
    for ii = 1:size(TAC, 2)
        for jj = 1:ii - 1
            concat = [TAC(:, ii) TACint(:, ii) -TAC(:, jj) -TACint(:, jj)];
            [~, ~, V] = svd(concat);
            CpInts = [CpInts abs(concat(:, 1:2) * V(1:2, end)) ...
                abs(-concat(:, 3:4) * V(3:4, end))];
        end
    end
    
    % Now choose CpInt that is closest to all others.
    differences = zeros(size(CpInts, 2), 1);
    for ii = 1:size(CpInts, 2)
        differences(ii) = norm(CpInts - CpInts(:, ii), 'fro');
    end
    [~, minInd] = min(differences); % all positives!
    bestCpInt = CpInts(:, minInd);
        
    
    
    
    
%     n = size(TAC, 2);
%     k = 1; %counter
%     Cpint = zeros(size(TAC, 1), size(TAC, 2), size(TAC, 2));
%     dim = n*(n-1)/2;
%     Cpintlist = zeros(18,dim); %same thing with Cpint but in 2-dim
%     %va = zeros(4,dim);      %the [a1i,a2i,a1j,a2j] in Fields
% 
%     for i = 1:n
%         cti = TAC(:,i);
%         intcti = TACint(:, i);
% 
%         %cpi = pdata(:,2);
%         %cpi = cumtrapz(pdata(:,1),pdata(:,2));
%         for j = i + 1:n
%             ctj = TAC(:, j);
%             intctj = TACint(:, j);
%             % find the right singular vector va corresponding to smallest
%             % right singular value of c4
%             c4 = [cti, intcti, -ctj, -intctj];
%             [~, ~, V] = svd(c4);
%             %va(:,k) = V(:,end);
% 
%             %evaluate the first guesses of Cpint
%             Cpint(:,i,j) = abs(V(1,end))*cti +abs(V(2,end))*intcti;
%             Cpintlist(:,k) = Cpint(:,i,j);
%             k = k+1;
%         end
%     end
% 
%     %find one cpint which minimize the summation of distances
%     dist = zeros(1,dim);
% 
%     for l = 1:dim
%         for m = 1:dim
%             disttemp = norm(Cpintlist(:,l)-Cpintlist(:,m));
%             dist(l) = dist(l) + disttemp;
%         end
%     end
% 
%     [~, minl] = min(dist);
%     bestCpInt = Cpintlist(:,minl);
end

