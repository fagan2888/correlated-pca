%%%Wrapper to generate the plots for the analysis of SVD algorithm under
%%%non-isotropic and also data-dependent nosie assumptions.


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

AlRange = unique(ceil(linspace(1 ,  3 * n, 100)));

FinalSubspaceError = zeros(num_trials, length(AlRange));
EstimatedSubspaces = cell(num_trials, length(AlRange));


SE_theory = zeros(size(AlRange));
SE_theory_temp = zeros(size(AlRange));


for mc = 1 : num_trials
    
    fprintf('MC iteration %d..\n', mc);
    
    %%Generate true data
    A = zeros(r, t_max);
    for jj = 1 : r
        A(jj, :) = -BoundL(jj) + ...
            2 * BoundL(jj) * rand(1, t_max);
    end
    L = P * A;
    
    %%Generate noise -- independent, non-isotropic, bounded
    %V = zeros(n, t_max);
    
    C = zeros(r, t_max);
    for jj = 1 : r
        C(jj, :) = -diag_entries_noise(jj) + ...
            2 * diag_entries_noise(jj) * rand(1, t_max);
    end
    
    V = B * C;
    
    %Generate data-dependent noise
    %b0 = 0.25;
%     beta = 100;
%     I = eye(n);
%     s = 0.05 * n;
%     rho = 1;
%     q = 1e-3;
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
%     
%     Y = L + W + V;
    
    
    %% Perform SVD for different values of \alpha and check accuracy
    %%parallelized to increase speed.
    
    for ii = 1 : length(AlRange)
        alpha = AlRange(ii);
        
        %Generate data-dependent noise
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
        
        EmpiricalCovariance = 1 / alpha * Y(:, 1: alpha) * Y(:, 1: alpha)';
        
        EstimatedSubspaces{mc, ii} = simpleEVD(EmpiricalCovariance, r);
        
        FinalSubspaceError(mc, ii) ...
            = Calc_SubspaceError(EstimatedSubspaces{mc, ii}, P);
        
        %%compute theoretical bounds
        %uncorrelated bounds
        Sigma_v = B * diag(flip(diag_entries_noise.^2 / 6)) * B';
        %Sigma_v = diag(flip(diag_entries_noise.^2 / 6));
        XX = P' * Sigma_v * P;
        lambda_vp_minus = min(eig(XX));
        YY = Sigma_v - P * XX * P';
        lambda_vrest_plus = max(eig(YY));
        lambda_p_pperp = norm(P_perp' * Sigma_v * P);
        lambda_v_plus = norm(Sigma_v);
        
        lambda_minus = min(BoundL)^2 / 6;
        lambda_plus = max(BoundL)^2 / 6;
        
        f = lambda_plus / lambda_minus;
        
        g = max(lambda_v_plus / lambda_minus, sqrt(lambda_v_plus / lambda_minus * f));
        
        d_uncor_alpha = g * sqrt(n * log(n)/ alpha);
        d_uncor_denom_alpha = f * sqrt((r + log(n)) / alpha);
        
        %         SE_theory(ii) = ((lambda_p_pperp / lambda_minus) + d_uncor_alpha) / ...
        %             (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
        %             d_uncor_alpha - d_uncor_denom_alpha);
        
        %correlated bounds
        %b_0 = beta / alpha;
        d_cor_alpha = q * f * sqrt(r * log(n) / alpha);
        d_cor_denom_alpha = f * sqrt((r + log(n)) / alpha);
        
        %         SE_theory(ii) = (3 * sqrt(b_0) * q * f + d_cor_alpha) / ...
        %             (1 - 3 * sqrt(b_0) * q * f -d_cor_alpha - d_cor_denom_alpha);
        
        d_alpha(ii) =  max([q * f * sqrt(r * log(n) / alpha), ...
            sqrt(lambda_v_plus / lambda_minus * f) * sqrt(r * log(n) / alpha), ...
            lambda_v_plus/ lambda_minus  * sqrt(r * log(n) / alpha)]);
        d_denom_alpha(ii) = f * sqrt((r + log(n)) / alpha);
        
        SE_theory(ii) = (lambda_p_pperp / lambda_minus + sqrt(b_0) * (2*q + q^2) * f + d_alpha(ii)) / ...
            (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
            sqrt(b_0) * (2*q + q^2) * f - d_alpha(ii) - d_denom_alpha(ii));
        
        SE_theory_temp(ii) = (lambda_p_pperp / lambda_minus + sqrt(b_0) * (2*q + q^2) * f) / ...
            (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
            sqrt(b_0) * (2*q + q^2) * f);
    end
end

%% Visulize results

% AlRange_r5 = AlRange;
% SE_theory_r5 = SE_theory;
% FinalSubspaceError_r5 = FinalSubspaceError;
% save('r_5.mat', 'AlRange_r5', 'SE_theory_r5', 'FinalSubspaceError_r5')

figure
%subplot(211)
plot(AlRange, mean(FinalSubspaceError, 1), 'b-');
hold
plot(AlRange, max(FinalSubspaceError, [], 1), 'r-');
% plot(AlRange, SE_theory, 'g*-')
%plot(AlRange, SE_theory_temp, 'ks-.')
axis tight
xlabel('alpha');
ylabel('SE');
legend('mean SE', 'max SE', 'Predicted SE bound');
title('n = 100, r = r_v = 10')
%title('y = l + w + v')
% subplot(212)
% plot(AlRange, d_alpha, 'bo-')
% hold
% plot(AlRange, d_denom_alpha, 'rs-')
% axis tight
% xlabel('alpha');
% ylabel('d');
% legend('d(alpha)', 'd_{denom}(alpha)');

toc