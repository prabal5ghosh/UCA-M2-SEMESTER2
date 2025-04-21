% mranc_synthetic_signals.m
% Multi-Reference Adaptive Noise Cancelling (MRANC)
% Section 2.1 Synthetic Signals: Effects of μ and M

clc; clear; close all;

% ----------------------
% Parameters
% ----------------------
N        = 1000;          % Number of samples
fs       = 500;           % Sampling frequency (Hz)
t        = (0:N-1)/fs;    % Time vector

% Desired signal: 5 Hz sinusoid with random phase
desired_freq    = 5;                  % Hz
phase_desired   = rand*2*pi;          % random phase
d = sin(2*pi*desired_freq*t + phase_desired);

% Interference: 50 Hz sinusoid with controlled input SIR
interf_freq     = 50;        % Hz
phase_interf    = rand*2*pi; % random phase
SIR_in          = 0;         % Input SIR (dB)
A_interf        = 10^(-SIR_in/20);
interference    = A_interf * sin(2*pi*interf_freq*t + phase_interf);

% Primary input: desired + interference
x = d + interference;

% Reference signals: pure interference components
z1 = sin(2*pi*interf_freq*t);
z2 = cos(2*pi*interf_freq*t);
Z  = [z1; z2];      % 2 × N matrix of references

% ----------------------
% Parameter Combinations
% ----------------------
mu_values = [0.0005, 0.005, 0.01, 0.01, 0.01];
M_values  = [0,       0,     5,    50,   550];

% Preallocate results: [mu, M, SIR_out, SIR_imp, conv_time]
results = nan(numel(mu_values), 5);

% ----------------------
% Loop Over Combinations
% ----------------------
for i = 1:numel(mu_values)
    mu = mu_values(i);
    M  = M_values(i);

    % Run MRANC
    [s_out, W] = mranc(x, Z, M, mu);

    % Steady-state start index
    conv_start = round(0.7*N) + 1;

    % Compute output SIR
    noise_residual = s_out(conv_start:end) - d(conv_start:end);
    SIR_out = 10*log10(var(d(conv_start:end)) / var(noise_residual));

    % Compute input SIR (actual)
    noise_orig = x(conv_start:end) - d(conv_start:end);
    SIR_in_actual = 10*log10(var(d(conv_start:end)) / var(noise_orig));

    % SIR improvement
    SIR_imp = SIR_out - SIR_in_actual;

    % Estimate convergence time
    err    = abs(s_out - d);
    err_s  = movmean(err, 50);
    thresh = 1.5 * mean(err_s(conv_start:end));
    idx    = find(err_s < thresh, 1);
    if isempty(idx)
        conv_time = NaN;
    else
        conv_time = idx;
    end

    % Store results
    results(i,:) = [mu, M, SIR_out, SIR_imp, conv_time];

    % ----------------------
    % Plotting
    % ----------------------
    figure;
    subplot(3,1,1);
    plot(t, d, 'b', t, x, 'r', t, s_out, 'g', 'LineWidth', 1.2);
    title(sprintf('MRANC Output (\\mu=%.4f, M=%d)', mu, M));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Desired', 'Input', 'Output'); grid on;

    subplot(3,1,2);
    plot(t, err, 'r', t, err_s, 'b', 'LineWidth', 1.2);
    xlabel('Time (s)'); ylabel('Error'); title('Error Signal'); grid on;
    if ~isnan(conv_time)
        hold on;
        plot([t(conv_time) t(conv_time)], ylim, 'k--', 'LineWidth', 1);
        text(t(conv_time), max(err_s), 'Convergence');
    end

    subplot(3,1,3);
    [Px, f] = pwelch(x, hamming(256), 128, 1024, fs);
    [Ps, ~] = pwelch(s_out, hamming(256), 128, 1024, fs);
    plot(f, 10*log10(Px), 'r', f, 10*log10(Ps), 'g', 'LineWidth', 1.2);
    xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title('Power Spectral Density');
    legend('Input', 'Output'); grid on;
    xlim([0 100]);
end

% ----------------------
% Display Results Table
% ----------------------
fprintf('-----------------------------------------------------------------\n');
fprintf('|   mu    |  M  |  SIR_out (dB)  |  SIR_imp (dB)  | Convergence |\n');
fprintf('-----------------------------------------------------------------\n');
for i = 1:size(results,1)
    fprintf('| %6.4f | %4d | %14.2f | %14.2f | %11.0f |\n', ...
        results(i,1), results(i,2), results(i,3), results(i,4), results(i,5));
end
fprintf('-----------------------------------------------------------------\n');

% ----------------------
% MRANC Function Definition
% ----------------------
function [s, W] = mranc(x, Z, M, mu)
    [R, N] = size(Z);
    W = zeros(R*(M+1), 1);
    s = zeros(1, N);

    % Build Toeplitz matrix
    ZZ = [];
    for r = 1:R
        ZZ = [ZZ; toeplitz([Z(r,1), zeros(1,M)], Z(r,:) )];
    end

    % Adaptive loop
    for n = M+1:N
        z_vec = ZZ(:, n);
        y     = W' * z_vec;
        s(n)  = x(n) - y;
        W     = W + 2*mu * s(n) * z_vec;
    end
end
