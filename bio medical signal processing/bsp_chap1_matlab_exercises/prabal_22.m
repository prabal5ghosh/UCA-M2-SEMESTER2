% fetal_ecg_extraction.m
% Adaptive cancellation of MECG to extract FECG from abdominal leads

clc; clear; close all;

% ----------------------
% Load ECG data
% ----------------------
% Replace the path below with your actual file path
data = load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\ecg_mother.mat');

% Inspect variable name if unsure:
disp('Variables in MAT file:');
disp(fieldnames(data));

% Use the correct variable name (update if different)
ECG = data.Y;  % Change to match your loaded variable name

% Normalize ECG signals (remove DC offset)
ECG = ECG - mean(ECG);

% ----------------------
% Parameters
% ----------------------
fs = 1000;        % Sampling frequency (Hz)
mu = 0.01;        % Step size for adaptation
M = 0;            % No time dispersion (simple LMS, single-tap)

% ----------------------
% Separate abdominal and thoracic signals
% ----------------------
abd = ECG(:,1:5);    % Abdominal (mixed FECG + MECG)
thr = ECG(:,6:8);    % Thoracic (MECG only)

[N, ~] = size(ECG);
fetal_estimates = zeros(N, 5);  % Output matrix for estimated FECG signals

% ----------------------
% MRANC function (no time dispersion, M = 0)
% ----------------------
function s = mranc_basic(x, Z, mu)
    [N, R] = size(Z);     % Z: Nx3 matrix (3 thoracic leads)
    W = zeros(R, 1);      % Adaptive weights
    s = zeros(N, 1);      % Output signal (filtered)
    for n = 1:N
        zvec = Z(n, :)';  % Current reference vector (3x1)
        y = W' * zvec;    % Estimated MECG
        s(n) = x(n) - y;  % Error = Abdominal - Estimated MECG
        W = W + 2 * mu * s(n) * zvec;  % Update weights
    end
end

% ----------------------
% Run adaptive cancellation on each abdominal lead
% ----------------------
for i = 1:5
    x = abd(:,i);           % One abdominal signal
    Z = thr;                % All 3 thoracic reference signals
    fetal_estimates(:,i) = mranc_basic(x, Z, mu);  % Run MRANC
end

% ----------------------
% Plot results for each abdominal channel
% ----------------------
time = (0:N-1)/fs;

for i = 1:5
    figure;
    subplot(3,1,1);
    plot(time, abd(:,i), 'r');
    title(sprintf('Abdominal Lead %d (Mixed FECG + MECG)', i));
    ylabel('Amplitude'); grid on;

    subplot(3,1,2);
    plot(time, thr(:,1), 'k'); % Plot 1 thoracic channel (as example)
    title('Thoracic Lead (MECG Reference)');
    ylabel('Amplitude'); grid on;

    subplot(3,1,3);
    plot(time, fetal_estimates(:,i), 'b');
    title(sprintf('Estimated FECG (Lead %d)', i));
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
end





%%

% ----------------------
% Analysis of Signal Estimation Quality
% ----------------------

% 2(a) Time-domain zoom (first 2 seconds)
t_zoom = time <= 2;  % Logical index for first 2 seconds

for i = 1:5
    figure;
    plot(time(t_zoom), abd(t_zoom,i), 'r', 'DisplayName', 'Abdominal Input');
    hold on;
    plot(time(t_zoom), fetal_estimates(t_zoom,i), 'b', 'DisplayName', 'Estimated FECG');
    title(sprintf('Lead %d: Convergence in First 2 Seconds', i));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend; grid on;
end

% 2(b) Frequency-domain analysis (Welch PSD from 0–50 Hz)
% We'll consider the steady-state part: last 30% of samples
ss_start = round(0.7 * N);
win = hamming(1024);
nfft = 1024;

for i = 1:5
    [Pxx_in, f]  = pwelch(abd(ss_start:end, i), win, [], nfft, fs);
    [Pxx_out, ~] = pwelch(fetal_estimates(ss_start:end, i), win, [], nfft, fs);

    % Limit to 0–50 Hz
    f_range = f <= 50;

    figure;
    plot(f(f_range), 10*log10(Pxx_in(f_range)), 'r', 'DisplayName', 'Input (Abdominal)');
    hold on;
    plot(f(f_range), 10*log10(Pxx_out(f_range)), 'b', 'DisplayName', 'Output (Estimated FECG)');
    title(sprintf('Lead %d: Welch Spectrum (0–50 Hz)', i));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend; grid on;
end
