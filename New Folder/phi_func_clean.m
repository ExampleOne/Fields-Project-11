
function phi_SA = phi_func_clean(LR, t, TACs, frameWidths)
% PHI_SA is the error function for the SIME simulated annealing
%   t should be a column vector.
    [Y1, names] = readTAC('/home/jouyang/Downloads/GIT_code/New Folder/fullTAC0sigma.tac');
    a = [702.9041; 12.8231; 9.2728];
    b = [-0.0182; -0.004; -0.0001];
    %upp = 5500; % upper bound
    mid = 136; % intermediate bound
    low = 76; % lower bound
    sample = 136; % one blood sampling time
    
    L = reshape(LR(:, 1), [2 5])';
    R = reshape(LR(:, 2), [2 5])';
    f_theta = zeros(28, 5); % 28 frames, 5 ROIs
    
    % exp_parts ignore the terms when x is the upper bound of +inf
    % for ii = 1:size(L, 2)
    for ii = 1:size(L, 1) % L is 5 by 2
        flat_part = 725 * L(ii, :) .* ...
            ( exp( - R(ii, :) .* (t - low)) -  exp( - R(ii, :) .* (t - mid) ) ) ...
            ./ R(ii, :);
        exp_part1 = a(1) * L(ii, :) .* exp( -R(ii, :) .* (t - mid)) ./ ...
            (b(1) - R(ii, :));
        exp_part2 = a(2) * L(ii, :) .* exp( -R(ii, :) .* (t - mid)) ./ ...
            (b(2) - R(ii, :));
        exp_part3 = a(3) * L(ii, :) .* exp( -R(ii, :) .* (t - mid)) ./ ...
            (b(3) - R(ii, :));
        conv = flat_part + exp_part1 + exp_part2 + exp_part3;
        f_theta(:, ii) = sum(conv, 2); % f_theta is 28 by 5
    end
    
    f_theta;
    figure
    plot(f_theta(:,1), t, Y1(:,3), t); 
    
    phi_SA = norm((f_theta - TACs) .* sqrt(frameWidths), 'fro') ^ 2;
end

% Draft:
    %   for j = 1:28
       
    %   For 1 regeion (L(1), L(2), R(1), R(2)):
    %   f_theta(j) = 725*LR(1,1)*(exp(LR(1,2)*(m-l))-exp(LR(1,2)*(t(j)-m)))/LR(1,2) ...
    %              + 725*LR(2,1)*(exp(LR(2,2)*(t(j)-l))-exp(LR(2,2)*(t(j)-m)))/LR(2,2) ...
    %              + sum(-LR(1,1).* a*exp(LR(1,2)*(t(j)-m))./(b-LR(1,2))) + sum(LR(1,1).* a*exp(LR(1,2)*(t(j)-u))./(b-LR(1,2))) ...
    %              + sum(-LR(2,1).* a*exp(LR(2,2)*(t(j)-m))./(b-LR(2,2))) + sum(LR(2,1).* a*exp(LR(2,2)*(t(j)-u)+b(3)*(u-s))./(b-LR(2,2)));

    %   For 5 regions (10 L's & 10 R's):
    %   f_theta(j) = 725 .* LR(1:2:9,1) .* (exp(LR(1:2:9,2)*(mid-low)) - exp(LR(1:2:9,2) *(t(j)-mid))) ./ LR(1:2:9,2) ...
    %              + 725 .* LR(2:2:10,1) .* (exp(LR(2:2:10,2) .* (t(j)-low)) - exp(LR(2:2:10,2) .* (t(j)-mid))) ./ LR(2:2:10,2) ...
    %              + sum(-LR(1:2:9,1) .* a .* exp(LR(1:2:9,2) .* (t(j)-mid))./(b - LR(1:2:9,2))) + sum(LR(1:2:9,1) .* a .* exp(LR(1:2:9,2) .* (t(j)-upp))./(b - LR(1:2:9,2))) ...
    %              + sum(-LR(2:2:10,1) .* a .* exp(LR(2:2:10,2) .* (t(j)-mid))./(b -LR(2:2:10,2))) + sum(LR(2:2:10,1) .* a .* exp(LR(2:2:10,2) .* (t(j)-upp)+b(3)*(upp-sample)) ./ (b-LR(2:2:10,2)));

    %         phi_SA = sum((w(j) * (Y1_Cerebellum(j) - f_theta(j))^2)) +  100 * (Pl-725)^2;
    %     end
