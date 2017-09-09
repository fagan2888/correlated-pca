%%%Wrapper to generate the plots for the analysis of SVD algorithm under
%%%non-isotropic and also data-dependent nosie assumptions.


clear
clc
close all

%% Initialization
tic
%n = 100;
%r = 1;

num_trials = 50;
alpha_num = 10;

% SE_theory = zeros(210, alpha_num);
% SE_theory_temp = zeros(210, alpha_num);
%
% mean_SE = zeros(210, alpha_num);
% max_SE = zeros(210, alpha_num);

SE_theory = zeros(10, alpha_num);
SE_theory_temp = zeros(10, alpha_num);

mean_SE = zeros(10, alpha_num);
max_SE = zeros(10, alpha_num);

all_errors = cell(10, 1);

cnt = 1;

for nn = 100
% for nn = unique(ceil(linspace(50, 1000, 10)))
    n = nn;
        for rr = [1 : 1 : 10]
%     for rr = 1
        r = rr;
        
        fprintf('n = %d, \t r = %d, \n', n, r);
        
        t_max = 500;
        
        U = orth(randn(n, n));
        P = U(:, 1 : r);
        P_perp = U(:, r+1 : end);
        B = orth(randn(n, r));
        
        %         BoundL = linspace(6, 6, r);
        %         diag_entries_noise = linspace(1, 1.1, r);
        
        
        AlRange = unique(ceil(linspace(10,  500, alpha_num)));
        
        FinalSubspaceError = zeros(num_trials, length(AlRange));
        EstimatedSubspaces = cell(num_trials, length(AlRange));
        
        
        for mc = 1 : num_trials
            
            
            %% Perform SVD for different values of \alpha and check accuracy
            %%parallelized to increase speed.
            
            parfor ii = 1 : length(AlRange)
                
                BoundL = linspace(6, 6, r);
                diag_entries_noise = linspace(1, 1.1, r);
                
                
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
                
                alpha = AlRange(ii);
                fprintf('MC %d..\t alpha %d\n', mc, alpha);
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
                for jj = 1 : t_max
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
                
                %correlated bounds
                %                 d_cor_alpha = q * f * sqrt(r * log(n) / alpha);
                %                 d_cor_denom_alpha = f * sqrt((r + log(n)) / alpha);
                
                d_alpha(ii) =  max([q * f * sqrt(r * log(n) / alpha), ...
                    sqrt(lambda_v_plus / lambda_minus * f) * sqrt(r * log(n) / alpha), ...
                    lambda_v_plus/ lambda_minus  * sqrt(r * log(n) / alpha)]);
                d_denom_alpha(ii) = f * sqrt((r + log(n)) / alpha);
                
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
lambda_v_plus = 0.2017;
lambda_minus = 6;
f = 1;

PhaseTrans = zeros(10, alpha_num);
%thresh = 0.01;
thresh = 1.8 * (lambda_v_plus / lambda_minus + b_0 * (2 * q + q^2) * f);
for ii = 1 : 10
    temp = all_errors{ii};
    for jj = 1 : alpha_num
        temp1 = temp(:, jj);
        PhaseTrans(ii, jj) = length(find(temp1 <= thresh));
    end
end

figure
imagesc(AlRange, [1 : 10], PhaseTrans);
xlabel('\alpha')
ylabel('n')
colormap('gray')

save('data/phase_trans_vs_r_mc50.mat')



% figure
% plot(AlRange, mean_SE(1, :), 'ro-')
% hold
% plot(AlRange, mean_SE(2, :), 'rs-.')
% plot(AlRange, mean_SE(3, :), 'r*--')
% plot(AlRange, mean_SE(4, :), 'go-')
% plot(AlRange, mean_SE(5, :), 'gs-.')
% plot(AlRange, mean_SE(6, :), 'g*--')
% plot(AlRange, mean_SE(7, :), 'bo-')
% plot(AlRange, mean_SE(8, :), 'bs-.')
% plot(AlRange, mean_SE(9, :), 'b*--')
% axis tight
% xlabel('alpha')
% ylabel('SE')
% title('mean SE')
% legend('n = 100, r=1', 'n = 100, r=5', 'n = 100, r=10', ...
%     'n = 500, r=1', 'n = 500, r=5', 'n = 500, r=10', ...
%     'n = 1000, r=1', 'n = 1000, r=5', 'n = 1000, r=10')
%
% figure
% plot(AlRange, max_SE(1, :), 'ro-')
% hold
% plot(AlRange, max_SE(2, :), 'rs-.')
% plot(AlRange, max_SE(3, :), 'r*--')
% plot(AlRange, max_SE(4, :), 'go-')
% plot(AlRange, max_SE(5, :), 'gs-.')
% plot(AlRange, max_SE(6, :), 'g*--')
% plot(AlRange, max_SE(7, :), 'bo-')
% plot(AlRange, max_SE(8, :), 'bs-.')
% plot(AlRange, max_SE(9, :), 'b*--')
% axis tight
% xlabel('alpha')
% ylabel('SE')
% title('max SE')
% legend('n = 100, r=1', 'n = 100, r=5', 'n = 100, r=10', ...
%     'n = 500, r=1', 'n = 500, r=5', 'n = 500, r=10', ...
%     'n = 1000, r=1', 'n = 1000, r=5', 'n = 1000, r=10')



% AlRange_r5 = AlRange;
% SE_theory_r5 = SE_theory;
% FinalSubspaceError_r5 = FinalSubspaceError;
% save('r_5.mat', 'AlRange_r5', 'SE_theory_r5', 'FinalSubspaceError_r5')

% figure
% subplot(211)
% plot(AlRange, mean(FinalSubspaceError, 1), 'bo-');
% hold
% plot(AlRange, max(FinalSubspaceError, [], 1), 'rs-');
% plot(AlRange, SE_theory, 'g*-')
% plot(AlRange, SE_theory_temp, 'ks-.')
% axis tight
% xlabel('alpha');
% ylabel('SE');
% legend('mean SE', 'max SE', 'Predicted SE bound', 'SE bound with d=0');
% title('n = 1000, r = r_v = 10')
% %title('y = l + w + v')
% subplot(212)
% plot(AlRange, d_alpha, 'bo-')
% hold
% plot(AlRange, d_denom_alpha, 'rs-')
% axis tight
% xlabel('alpha');
% ylabel('d');
% legend('d(alpha)', 'd_{denom}(alpha)');

toc