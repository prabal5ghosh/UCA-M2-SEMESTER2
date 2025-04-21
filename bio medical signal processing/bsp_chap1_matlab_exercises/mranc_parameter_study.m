% mranc_parameter_study.m
% MRANC Parameter Study in MATLAB

clc; clear; close all;

% ----------------------
% Generate synthetic desired and interference signals
% ----------------------
N     = 1000;          % Number of samples
fs    = 500;           % Sampling frequency (Hz)
t     = (0:N-1)/fs;    % Time vector

% Desired signal: 5 Hz sinusoid with random phase
desired = sin(2*pi*5*t + rand*2*pi);

% Interference: 50 Hz sinusoid scaled to input SIR of 0 dB
SIR_in_dB = 0;
A_int     = 10^(-SIR_in_dB/20);
interf    = A_int * sin(2*pi*50*t + rand*2*pi);

% Primary input
x = desired + interf;

% Reference signals
z1 = sin(2*pi*50*t);
z2 = cos(2*pi*50*t);
Z  = [z1; z2];      % 2-by-N reference matrix

% ----------------------
% Parameter combinations
% ----------------------
mu_vals = [0.0005, 0.005, 0.01, 0.01, 0.01];
M_vals  = [0,       0,     5,    50,   550];
numCases = numel(mu_vals);

% Preallocate results
SIR_out      = nan(1,numCases);
SIR_imp      = nan(1,numCases);
conv_samples = nan(1,numCases);

% ----------------------
% MRANC function (nested)
% ----------------------
function s = mranc_algo(x, Z, M, mu)
    [R, N] = size(Z);
    % Build Toeplitz matrix of references
    ZZ = [];
    for r = 1:R
        ZZ = [ZZ; toeplitz([Z(r,1) zeros(1,M)], Z(r,:))];
    end
    W = zeros(R*(M+1),1);
    s = zeros(1,N);
    for n = M+1:N
        zvec = ZZ(:,n);
        y    = W' * zvec;
        s(n) = x(n) - y;
        W    = W + 2*mu * s(n) * zvec;
    end
end

% ----------------------
% Loop through parameter sets
% ----------------------
for k = 1:numCases
    mu = mu_vals(k);
    M  = M_vals(k);
    % Run MRANC
    s_out = mranc_algo(x, Z, M, mu);
    % Steady-state region
    ss_start = round(0.7 * N) + 1;
    % Compute output SIR
    res_noise = s_out(ss_start:end) - desired(ss_start:end);
    SIR_out(k) = 10*log10(var(desired(ss_start:end)) / var(res_noise));
    % Compute actual input SIR
    orig_noise = x(ss_start:end) - desired(ss_start:end);
    SIR_in_actual = 10*log10(var(desired(ss_start:end)) / var(orig_noise));
    SIR_imp(k) = SIR_out(k) - SIR_in_actual;
    % Estimate convergence time
    err      = abs(s_out - desired);
    err_s    = movmean(err, 50);
    thresh   = 1.5 * mean(err_s(ss_start:end));
    idx_conv = find(err_s < thresh, 1);
    if ~isempty(idx_conv)
        conv_samples(k) = idx_conv;
    end
end

% ----------------------
% Display results in table
% ----------------------
fprintf('-------------------------------------------------------------\n');
fprintf('|   mu    |   M   |  SIR_out (dB)  |  SIR_imp (dB)  | Conv samples |\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:numCases
    fprintf('| %7.4f | %5d | %14.2f | %14.2f | %12.0f |\n', ...
        mu_vals(k), M_vals(k), SIR_out(k), SIR_imp(k), conv_samples(k));
end
fprintf('-------------------------------------------------------------\n');

% End of mranc_parameter_study.m


% Choose a case to plot
k_plot = 3;
mu = mu_vals(k_plot);
M  = M_vals(k_plot);
s_out = mranc_algo(x, Z, M, mu);

figure;
plot(t, desired, 'g', 'LineWidth', 1.2); hold on;
plot(t, x, 'r', 'LineWidth', 1);
plot(t, s_out, 'b', 'LineWidth', 1);
legend('Desired', 'Noisy Input', 'Filtered Output');
xlabel('Time (s)');
ylabel('Amplitude');
title(sprintf('Signal Comparison (\\mu = %.4f, M = %d)', mu, M));
grid on;


figure;
plot(t, abs(s_out - desired), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Absolute Error');
title(sprintf('Error Signal |s(n) - d(n)| (\\mu = %.4f, M = %d)', mu, M));
grid on;


figure;
plot(M_vals, SIR_imp, 'o-b', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Filter Order (M)');
ylabel('SIR Improvement (dB)');
title('Effect of Filter Order on SIR Improvement');
grid on;




