% Brain regions (ROIs) to consider: Cerebellum, Temporal, Putamen,
% Thalamus, Cingulate (i = 5)

% concentration in region i:
% for ii = 1:14
Cp = dlmread('/home/jouyang/Downloads/GIT_code/New Folder/PlasmaGenerated.smpl', '\t', 1, 0);
[Y1, names] = readTAC('/home/jouyang/Downloads/GIT_code/New Folder/fullTAC0sigma.tac');
Cp_real = dlmread('/home/jouyang/Downloads/GIT_code/New Folder/pabloModel_0sigma.smpl', '\t', 1, 0);
% Cp_real = dlmread(['Users/shixiaoying/Desktop/MATLAB/newtacs/plasma' num2str(ii) '.smpl'], '\t', 1, 0);

TACs = Y1(:, [3 4 7 8 9]);
names = names([3 4 7 8 9]);
t = Y1(:,2);
Cp_real_t = Cp_real(:,1);
k = find(Cp_real_t < 136.5 & 135.5 < Cp_real_t);
Pl = Cp_real(k,2);
w = Y1(:, 2) - Y1(:, 1); % time frame w
% w(j) = (Y1(j, 2) - Y1(j, 1))';


% t_i = 0;
% t_f = Y1(29,2);
% t_p = Cp(:,1);

% Other attempts for doing convolution:
% f = @(t, R, L) simple_conv(L(1) * exp(-R(1) * t),Cp) + simple_conv(L(2)*exp(-R(2) * t), Cp);
% f = @(t, R, L) conv(L(1) * exp(-R(1) * t),Cp) + conv(L(2)*exp(-R(2) * t), Cp); % K = 2

% Skip:
% To estimate ROI-specific parameters (i.e. L(i), R(i))
% minimization using iterative nonlinear least square algorithm over all Li and Ri
% s = lsqnonlin(sumwYf, t)


% Simultaneous Estimation:
% minimization using simulated annealing that outputs parameters L(ik), R(ik), and f_theta for ROI i=1:5 and
% tissue compartment k=1:2

% minimizing objective function phi:
tic;
LR = simulatedAnnealing(@(LR) phi_func_clean(LR, t, TACs, w), ...
    [1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001; 1, 0.001], ...
    [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0], ...
    [+1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2; ...
    +1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2; +1e1, +1e-2], 1e10, 1e-5, 0.8);
toc;
display(LR);
display(phi_func_clean(LR, t, TACs, w));

% Vt for Cerebellum from this TAC1 
V_t_1 = sum(LR(1:2,1) ./ LR(1:2,2));
V_t_2 = sum(LR(3:4,1) ./ LR(3:4,2));
V_t_3 = sum(LR(5:6,1) ./ LR(5:6,2));
V_t_4 = sum(LR(7:8,1) ./ LR(7:8,2));
V_t_5 = sum(LR(9:10,1) ./ LR(9:10,2));

disp(['Vt of Cerebellar cortex = ' num2str(V_t_1)]);
disp(['Vt of Temporal cortex = ', num2str(V_t_2)]);
disp(['Vt of Putamen = ', num2str(V_t_3)]);
disp(['Vt of Thalamus = ', num2str(V_t_4)]);
disp(['Vt of Anterior = ', num2str(V_t_5)]);

% disp(['Vt of Cerebellum for plasma' num2str(ii) ' : ' ]);
    % display(V_t(:, ii)');
% end


% Future work:
% (sort of done) Calculate better bounds for L's and R's on paper, by pen
% (done) Split off convolution (f_theta) into another file
% Introduce a piece-wise function (Cp_theta) to estimate
% Add graphs of sourceCp (Cp_real) v.s. generatedCp (from SIME); inputTACs v.s. generatedTACs(from f_theta) 
% Use TACs with noise instead of fullTAC0sigma and judge the plots
% Use patients' TACs and determine the performance based on 2 plots



% Reference:
% constr is an implementation of sequential quadratic programming nonlinear minimization algorithm in Matlab 5.2.    
% model3Exp = @(c,x) c(1)*exp(-c(2)*x) .* x + c(3)*exp(-c(4)*x) +c(5)*exp(-(x-c(6)) .^ 2 * c(7));
% phi3Exp = @(c) norm(model3Exp(c, data(:, 1)) - data(:, 2)); %plasmaPhi(data, model, c);





