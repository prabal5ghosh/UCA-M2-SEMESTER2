% mranc_leakage_effects.m
% Study the effect of desired signal leakage into reference signals

clc; clear; close all;

% Generate synthetic desired and interference signals

N = 1000;
fs = 500;
t = (0:N-1)/fs;

desired = sin(2*pi*5*t + rand*2*pi); % Desired 5Hz
SIR_in_dB = 0;
A_int = 10^(-SIR_in_dB/20);
interf = A_int * sin(2*pi*50*t + rand*2*pi); % Interference 50Hz

x = desired + interf; % Primary input

% Contaminated reference signals (10 dB attenuation of desired)
atten_dB = 10;
alpha = 10^(-atten_dB/20); % Attenuation factor
z1 = sin(2*pi*50*t) + alpha * desired;
z2 = cos(2*pi*50*t) + alpha * desired;
Z = [z1; z2]; % Contaminated reference matrix

% ----------------------
% Parameter combinations
% ----------------------
mu_vals = [0.0005, 0.005, 0.01, 0.01, 0.01];
M_vals  = [0,       0,     5,    50,   550];
numCases = numel(mu_vals);

% Preallocate
SIR_out = nan(1,numCases);
SIR_imp = nan(1,numCases);
conv_samples = nan(1,numCases);

% MRANC function
function s = mranc_algo(x, Z, M, mu)
    [R, N] = size(Z);
    ZZ = [];
    for r = 1:R
        ZZ = [ZZ; toeplitz([Z(r,1) zeros(1,M)], Z(r,:))];
    end
    W = zeros(R*(M+1),1);
    s = zeros(1,N);
    for n = M+1:N
        zvec = ZZ(:,n);
        y = W' * zvec;
        s(n) = x(n) - y;
        W = W + 2 * mu * s(n) * zvec;
    end
end

% Run experiments for each parameter set
for k = 1:numCases
    mu = mu_vals(k);
    M = M_vals(k);
    s_out = mranc_algo(x, Z, M, mu);

    % Steady state region
    ss_start = round(0.7 * N) + 1;
    res_noise = s_out(ss_start:end) - desired(ss_start:end);
    SIR_out(k) = 10*log10(var(desired(ss_start:end)) / var(res_noise));
    orig_noise = x(ss_start:end) - desired(ss_start:end);
    SIR_in_actual = 10*log10(var(desired(ss_start:end)) / var(orig_noise));
    SIR_imp(k) = SIR_out(k) - SIR_in_actual;

    % Convergence time estimate
    err = abs(s_out - desired);
    err_s = movmean(err, 50);
    thresh = 1.5 * mean(err_s(ss_start:end));
    idx_conv = find(err_s < thresh, 1);
    if ~isempty(idx_conv)
        conv_samples(k) = idx_conv;
    end
end

% Display results
fprintf('\nEffects of Reference Contamination (10 dB Desired Signal Leakage):\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('|   mu    |   M   |  SIR_out (dB)  |  SIR_imp (dB)  | Conv samples |\n');
fprintf('--------------------------------------------------------------------------\n');
for k = 1:numCases
    fprintf('| %7.4f | %5d | %14.2f | %14.2f | %12.0f |\n', ...
        mu_vals(k), M_vals(k), SIR_out(k), SIR_imp(k), conv_samples(k));
end
fprintf('--------------------------------------------------------------------------\n');

% Plotting
% ----------------------
colors = lines(numCases);
figure;
for k = 1:numCases
    mu = mu_vals(k);
    M = M_vals(k);
    s_out = mranc_algo(x, Z, M, mu);
    subplot(numCases,1,k);
    plot(t, desired, 'g', 'LineWidth', 1); hold on;
    plot(t, x, 'r:', 'LineWidth', 1);
    plot(t, s_out, 'b', 'LineWidth', 1);
    title(sprintf('\\mu = %.4f, M = %d', mu, M));
    ylabel('Amplitude');
    if k == 1
        legend('Desired', 'Input', 'Output');
    end
end
xlabel('Time (s)');
sgtitle('MRANC with Contaminated Reference Signals (Leakage 10 dB)');

% Error plots
figure;
for k = 1:numCases
    mu = mu_vals(k);
    M = M_vals(k);
    s_out = mranc_algo(x, Z, M, mu);
    err = abs(s_out - desired);
    subplot(numCases,1,k);
    plot(t, err, 'LineWidth', 1.2);
    ylabel('|e(n)|');
    title(sprintf('Error Signal \\mu = %.4f, M = %d', mu, M));
end
xlabel('Time (s)');
sgtitle('Error Signal |s(n) - d(n)|');

% SIR improvement vs M
figure;
plot(M_vals, SIR_imp, 'o-b', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Filter Order (M)');
ylabel('SIR Improvement (dB)');
title('Effect of Filter Order on SIR Improvement (Contaminated References)');
grid on;
