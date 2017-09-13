%%This file contains the demo to validate the bound for the SVD guarantees
%%as predicted by our analysis for the bounded noise case. The result shows
%%that the bound obtained by our theorem is a good approximation of the
%%true (mean and max) subspace errors numerically for a wide range of
%%alpha.


clear
clc
close all

%% Initialization
tic
n = 100;
t_max = 7000;
r = 1;
U = orth(randn(n, n));
P = U(:, 1 : r);
P_perp = U(:, r+1 : end);
B = orth(randn(n, r));

BoundL = linspace(6, 6, r);
diag_entries_noise = linspace(1, 1.1, r);
num_trials = 3;

AlRange = unique(ceil(linspace(50 ,  3 * n, 100))); 
%choose AlRange smartly, otherwise the bounds will be meaningless

FinalSubspaceError = zeros(num_trials, length(AlRange));
EstimatedSubspaces = cell(num_trials, length(AlRange));


SE_theory = zeros(size(AlRange));
SE_theory_temp = zeros(size(AlRange));


for mc = 1 : num_trials
    
    fprintf('MC iteration %d..\n', mc);

    %% Data Generation
    A = zeros(r, t_max);
    for jj = 1 : r
        A(jj, :) = -BoundL(jj) + ...
            2 * BoundL(jj) * rand(1, t_max);
    end
    L = P * A;
    
    %%Generate noise -- independent, non-isotropic, bounded
    C = zeros(r, t_max);
    for jj = 1 : r
        C(jj, :) = -diag_entries_noise(jj) + ...
            2 * diag_entries_noise(jj) * rand(1, t_max);
    end
    
    V = B * C;
    
    for ii = 1 : length(AlRange) %can be replaced with parfor; 
                                 %then, better if other data generation 
                                 %is done inside too
        alpha = AlRange(ii);
        
        %Generate data-dependent noise -- doing in inner loop because the
        %noise (outlier-row-frac) depends on alpha
        b_0 = 0.05;
        beta = ceil(b_0 * alpha);
        I = eye(n);
        s = 0.05 * n;
        rho = 1;
        q = 1e-3;
        
        num_changes = floor(t_max/beta);
        T = zeros(n, t_max);
        for jj = 1 : num_changes
            bind = max(mod(floor((jj-1) * s/rho + 1), n), 1);
            sind = min(bind - 1 + s, n);
            idx = bind : sind;
            T(idx, (jj-1) * beta + 1 : jj * beta) = 1;
        end
        
        W = zeros(n, t_max);
        for jj = 1 : alpha
            idx = find(T(:, jj));
            temp = abs(randn(length(idx), n));
            Mst = q * temp / norm(temp * P);
            W(:, jj) = I(:, idx) * (Mst * L(:, jj));
        end
        Y = L + W + V;
        
        %% Perform SVD for different values of \alpha and check accuracy
        EmpiricalCovariance = 1 / alpha * Y(:, 1: alpha) * Y(:, 1: alpha)';
        EstimatedSubspaces{mc, ii} = simpleEVD(EmpiricalCovariance, r);
        FinalSubspaceError(mc, ii) ...
            = Calc_SubspaceError(EstimatedSubspaces{mc, ii}, P);
        
        %% Compute theoretical bounds
        %uncorrelated bounds
        Sigma_v = B * diag(flip(diag_entries_noise.^2 / 6)) * B';
        XX = P' * Sigma_v * P;
        lambda_vp_minus = min(eig(XX));
        YY = Sigma_v - P * XX * P';
        lambda_vrest_plus = max(eig(YY));
        lambda_p_pperp = norm(P_perp' * Sigma_v * P);
        lambda_v_plus = norm(Sigma_v);
        
        lambda_minus = min(BoundL)^2 / 6;
        lambda_plus = max(BoundL)^2 / 6;
        
        f = lambda_plus / lambda_minus;
        
        g = max(lambda_v_plus / lambda_minus, ...
            sqrt(lambda_v_plus / lambda_minus * f));
        
        %correlated bounds
        d_cor_alpha = q * f * sqrt(r * log(n) / alpha);
        d_cor_denom_alpha = f * sqrt((r + log(n)) / alpha);
        
        %final bounds
        d_alpha(ii) =  max([q * f * sqrt(r * log(n) / alpha), ...
            sqrt(lambda_v_plus / lambda_minus * f) * sqrt(r * log(n) / alpha), ...
            lambda_v_plus/ lambda_minus  * sqrt(r * log(n) / alpha)]);
        d_denom_alpha(ii) = f * sqrt((r + log(n)) / alpha);
        
        SE_theory(ii) = (lambda_p_pperp / lambda_minus + sqrt(b_0) * (2*q + q^2) * f + d_alpha(ii)) / ...
            (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
            sqrt(b_0) * (2*q + q^2) * f - d_alpha(ii) - d_denom_alpha(ii));
    end
end

%% Visulize results
figure
plot(AlRange, mean(FinalSubspaceError, 1), 'bo-');
hold
plot(AlRange, max(FinalSubspaceError, [], 1), 'rs-');
plot(AlRange, SE_theory, 'g*-')
axis tight
xlabel('alpha');
ylabel('SE');
legend('mean SE', 'max SE', 'Predicted SE bound');
title('n = 100, r = r_v = 10')
toc