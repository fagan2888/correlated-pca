%%%Wrapper to generate the plots for the analysis of SVD algorithm under
%%%non-isotropic and also data-dependent nosie assumptions.


clear
clc
%close all

%% Initialization
tic
n = 100;
t_max = 10000;
r = 10;
U = orth(randn(n, n));
P = U(:, 1 : r);
P_perp = U(:, r+1 : end);

BoundL = 100;
diag_entries_noise = linspace(1, 1.5, n);
num_trials = 100;

AlRange = unique(ceil(linspace(n * log(n) , 10 * n, 80)));

FinalSubspaceError = zeros(num_trials, length(AlRange));
EstimatedSubspaces = cell(num_trials, length(AlRange));


SE_theory = zeros(size(AlRange));


for mc = 1 : num_trials
    
    %%Generate true data
    A = -BoundL + 2 * BoundL * rand(r, t_max);
    L = P * A;
    
    %%Generate noise -- independent, non-isotropic, bounded
    V = zeros(n, t_max);
    for jj = 1 : n
        V(jj, :) = -diag_entries_noise(jj) + ...
            2 * diag_entries_noise(jj) * rand(1, t_max);
    end
    
    %%Generate data-dependent noise
    %     b0 = 0.25;
    %     beta = ceil(b0 * 1000);
    %     I = eye(n);
    %     s = 0.05 * n;
    %     rho = 1;
    %     q = 0.01;
    %
    %     num_changes = floor(t_max/beta);
    %     T = zeros(n, t_max);
    %     for jj = 1 : num_changes
    %         bind = max(mod(floor((jj-1) * s/rho + 1), n), 1);
    %         sind = min(bind - 1 + s, n);
    %         idx = bind : sind;
    %         T(idx, (jj-1) * beta + 1 : jj * beta) = 1;
    %     end
    %
    %     for jj = 1 : t_max
    %         idx = find(T(:, jj));
    %         temp = abs(randn(length(idx), n));
    %         Mst = q * temp / norm(temp * P);
    %         W(:, jj) = I(:, idx) * (Mst * L(:, jj));
    %     end
    
    Y = L + V;
    
    
    %% Perform SVD for different values of \alpha and check accuracy
    %%parallelized to increase speed.
    
    parfor ii = 1 : length(AlRange)
        alpha = AlRange(ii);
        
        EmpiricalCovariance = 1 / alpha * Y(:, 1: alpha) * Y(:, 1: alpha)';
        
        EstimatedSubspaces{mc, ii} = simpleEVD(EmpiricalCovariance, r);
        
        FinalSubspaceError(mc, ii) ...
            = Calc_SubspaceError(EstimatedSubspaces{mc, ii}, P);
        
        %%compute theoretical bounds
        Sigma_v = diag(flip(diag_entries_noise.^2 / 6));
        XX = P' * Sigma_v * P;
        lambda_vp_minus = min(eig(XX));
        YY = Sigma_v - P * XX * P';
        lambda_vrest_plus = max(eig(YY));
        lambda_p_pperp = norm(P_perp' * Sigma_v * P);
        lambda_v_plus = norm(Sigma_v);
        
        lambda_minus = BoundL^2 / 6;
        lambda_plus = BoundL^2 / 6;
        
        f = lambda_plus / lambda_minus;
        
        g = max(lambda_v_plus / lambda_minus, sqrt(lambda_v_plus / lambda_minus * f));
        d_alpha = g * sqrt(n * log(n)/ alpha);
        d_denom_alpha = f * sqrt((r + log(n)) / alpha);
        
        SE_theory(ii) = ((lambda_p_pperp / lambda_minus) + d_alpha) / ...
            (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - d_alpha - d_denom_alpha);
    end
end

%% Visulize results
figure
plot(AlRange, mean(FinalSubspaceError, 1));
hold
plot(AlRange, max(FinalSubspaceError, [], 1), 'r');
plot(AlRange, SE_theory, 'g')
axis tight
toc
