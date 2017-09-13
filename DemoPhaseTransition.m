%%Demo to generate the phase transition matrix to observe variation with
%%either (i) rank or (ii) Signal dimension

clear
clc
close all

%% Initialization
tic

num_trials = 5;
alpha_num = 10;

% nrange = unique(ceil(linspace(50, 500, 10)));
rrange = [1 : 10];

nrange = 100;
%rrange = 10;

SE_theory = zeros(length(rrange), alpha_num); %change to nrange if needed
SE_theory_temp = zeros(length(rrange), alpha_num);

mean_SE = zeros(length(rrange), alpha_num);
max_SE = zeros(length(rrange), alpha_num);

all_errors = cell(length(rrange), 1);

cnt = 1;

for nn = 100
    % for nn = nrange
    n = nn;
    for rr = rrange
        %     for rr = 5
        r = rr;
        
        fprintf('n = %d, \t r = %d, \n', n, r);
        t_max = 500;
        
        U = orth(randn(n, n));
        P = U(:, 1 : r);
        P_perp = U(:, r+1 : end);
        B = orth(randn(n, r));
        
        AlRange = unique(ceil(linspace(50,  t_max, alpha_num))); %choose smartly
        
        FinalSubspaceError = zeros(num_trials, length(AlRange));
        EstimatedSubspaces = cell(num_trials, length(AlRange));
        
        
        for mc = 1 : num_trials
            %%parallelized to increase speed.
            parfor ii = 1 : length(AlRange) %%if using parfor, note the bound parameters
                
                alpha = AlRange(ii);
                BoundL = linspace(6, 6, r);
                diag_entries_noise = linspace(0.5, 0.9, r);
                
                %% Data Generation
                %%bounded
                A = zeros(r, t_max);
                for jj = 1 : r
                    A(jj, :) = -BoundL(jj) + ...
                        2 * BoundL(jj) * rand(1, t_max);
                end
                
                %%gaussian
                %                 A = zeros(r, t_max);
                %                 for jj = 1 : r
                %                     A(jj, :) = 6 * randn(1, t_max);
                %                 end
                
                L = P * A;
                
                
                %%Generate anisotropic noise
                %%bounded
                V = zeros(n, t_max);
                C = zeros(r, t_max);
                for jj = 1 : r
                    C(jj, :) = diag_entries_noise(jj) * rand(1, t_max);
                end
                V = B * C;
                
                %%gaussian
                %                 V = zeros(n, t_max);
                %                 C = zeros(r, t_max);
                %                 for jj = 1 : r
                %                     C(jj, :) = diag_entries_noise(jj) * rand(1, t_max);
                %                 end
                %                 V = B * C;
                
                fprintf('MC %d..\t alpha %d\n', mc, alpha);
                
                %Generate data-dependent noise
                %%%Generating support set and sparse vectors
                b_0 = 0.01;
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
                for jj = 1 : t_max
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
                
                %%compute theoretical bounds
                %uncorrelated bounds
                Sigma_v = B * diag(flip(diag_entries_noise.^2)) * B';
                XX = P' * Sigma_v * P;
                lambda_vp_minus = min(eig(XX));
                YY = Sigma_v - P * XX * P';
                lambda_vrest_plus = max(eig(YY));
                lambda_p_pperp = norm(P_perp' * Sigma_v * P);
                lambda_v_plus = norm(Sigma_v);
                
                %%bounded
                lambda_minus = min(BoundL)^2/3;
                lambda_plus = max(BoundL)^2/3;
                
                %%gaussian
                %                 lambda_minus = BoundL^2;
                %                 lambda_plus = BoundL^2;
                
                f = lambda_plus / lambda_minus;
                
                g = max(lambda_v_plus / lambda_minus, sqrt(lambda_v_plus / lambda_minus * f));
                
                %correlated bounds
                %                 d_cor_alpha = q * f * sqrt(r * log(n) / alpha);
                %                 d_cor_denom_alpha = f * sqrt((r + log(n)) / alpha);
                
                d_alpha(ii) =  max([q * f * sqrt(r * log(n) / alpha), ...
                    sqrt(lambda_v_plus / lambda_minus * f) * sqrt(r * log(n) / alpha), ...
                    lambda_v_plus/ lambda_minus  * sqrt(r * log(n) / alpha)]);
                
                d_denom_alpha(ii) = f * sqrt((r + log(n)) / alpha);
                
                %%gaussian
                %                 d_denom_alpha(ii) = g * sqrt(n / alpha);
                
                SE_theory(cnt, ii) = (lambda_p_pperp / lambda_minus + sqrt(b_0) * (2*q + q^2) * f + d_alpha(ii)) / ...
                    (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
                    sqrt(b_0) * (2*q + q^2) * f - d_alpha(ii) - d_denom_alpha(ii));
                
                SE_theory_temp(cnt, ii) = (lambda_p_pperp / lambda_minus + sqrt(b_0) * (2*q + q^2) * f) / ...
                    (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus - ...
                    sqrt(b_0) * (2*q + q^2) * f);
            end
            
        end
        all_errors{cnt} = FinalSubspaceError;
        mean_SE(cnt, :) = mean(FinalSubspaceError, 1);
        max_SE(cnt, :) = max(FinalSubspaceError, [], 1);
        cnt = cnt + 1;
    end
end
%% Visulize results
b_0 = 0.05;
q = 0.001;
lambda_v_plus = 0.8100;
lambda_minus = 12;
f = 1;

% Sigma_v = B * diag(flip(diag_entries_noise.^2)) * B';
% %Sigma_v = diag(flip(diag_entries_noise.^2 / 6));
% XX = P' * Sigma_v * P;
% lambda_vp_minus = min(eig(XX));
% YY = Sigma_v - P * XX * P';
% lambda_vrest_plus = max(eig(YY));
% lambda_p_pperp = norm(P_perp' * Sigma_v * P);
% lambda_v_plus = norm(Sigma_v);
%
% thresh = 1.1 * (3 * sqrt(b_0) * q * f + ...
%     (lambda_p_pperp / lambda_minus) / (1 - (lambda_vrest_plus - lambda_vp_minus) / lambda_minus));

PhaseTrans = zeros(length(rrange), alpha_num);

thresh = 0.07; %%just set temporarily -- can adjust
for ii = 1 : length(rrange)
    temp = all_errors{ii};
    for jj = 1 : alpha_num
        temp1 = temp(:, jj);
        PhaseTrans(ii, jj) = length(find(temp1 <= thresh));
    end
end

figure
imagesc(AlRange, rrange, PhaseTrans);
xlabel('\alpha')
ylabel('n')
colormap('gray')

toc