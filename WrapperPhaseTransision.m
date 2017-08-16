%%% Code to generate the phase transition plot for PCA when the noise is
%%% non-isotropic.

clear
clc

%% Generate data
tic
n = 100;
t_max = 1000;
r = 5;
P = orth(randn(n, r));
BoundL = 5;
diag_entries_noise = linspace(1, 5, r);
num_trials = 100;

AlRange = unique(ceil(linspace(1, 10 * n, 80)));

FinalSubspaceError = zeros(num_trials, length(AlRange));
EstimatedSubspaces = cell(num_trials, length(AlRange));



for mc = 1 : num_trials
    A = -BoundL + 2 * BoundL * rand(r, t_max);
    L = P * A;
    
    %BoundV = BoundL;
    
    V = zeros(n, t_max);
    
    for ii = 1 : r
        V(ii, :) = -diag_entries_noise(ii) + ...
            2 * diag_entries_noise(ii) * rand(1, t_max);
    end
    
    Y = L + V;
    
    
    %% Perform SVD for different values of \alpha and check accuracy
    
    %AlRange = 50 : 50 : 1000;
    
    
    parfor ii = 1 : length(AlRange)
        alpha = AlRange(ii);
        
        EmpiricalCovariance = 1 / alpha * Y(:, 1: alpha) * Y(:, 1: alpha)';
        
        EstimatedSubspaces{mc, ii} = simpleEVD(EmpiricalCovariance, r);
        
        FinalSubspaceError(mc, ii) ...
            = Calc_SubspaceError(EstimatedSubspaces{mc, ii}, P);
    end
end

%% Visulize results
figure
plot(AlRange, mean(FinalSubspaceError, 1));

toc


