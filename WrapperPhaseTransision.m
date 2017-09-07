%%%Wrapper to generate the plots for the analysis of SVD algorithm under
%%%non-isotropic and also data-dependent nosie assumptions.


clear
clc
close all

%% Initialization
tic
n = 100;
t_max = 10000;
r = 30;
P = orth(randn(n, r));
BoundL = 5;
diag_entries_noise = linspace(1, 5, r);
num_trials = 100;
% W = zeros(n, t_max);
% V = zeros(n, t_max);

AlRange = unique(ceil(linspace(1, 10 * n, 80)));

FinalSubspaceError = zeros(num_trials, length(AlRange));
EstimatedSubspaces = cell(num_trials, length(AlRange));



for mc = 1 : num_trials
    
    %%Generate true data
    A = -BoundL + 2 * BoundL * rand(r, t_max);
    L = P * A;
    
    %%Generate noise -- independent, but anisotropic
    V = zeros(n, t_max);
    for jj = 1 : r
        V(jj, :) = -diag_entries_noise(jj) + ...
            2 * diag_entries_noise(jj) * rand(1, t_max);
    end
    
    %%Generate data-dependent noise
    b0 = 0.25;
    beta = ceil(b0 * 1000);
    I = eye(n);
    s = 0.05 * n;
    rho = 1;
    q = 0.01; 
    
    num_changes = floor(t_max/beta);
    T = zeros(n, t_max);
    for jj = 1 : num_changes
        bind = max(mod(floor((jj-1) * s/rho + 1), n), 1);
        sind = min(bind - 1 + s, n);
        idx = bind : sind;
        T(idx, (jj-1) * beta + 1 : jj * beta) = 1;
    end
    
    for jj = 1 : t_max
        idx = find(T(:, jj));
        temp = abs(randn(length(idx), n));
        Mst = q * temp / norm(temp * P);
        W(:, jj) = I(:, idx) * (Mst * L(:, jj));
    end
    
    Y = L + W + V;
    
    
    %% Perform SVD for different values of \alpha and check accuracy
    %%parallelized to increase speed.
    
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
hold
plot(AlRange, max(FinalSubspaceError, [], 1), 'r');
toc
