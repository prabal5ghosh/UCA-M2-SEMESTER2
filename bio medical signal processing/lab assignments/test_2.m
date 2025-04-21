% Load the data
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\ecg_mother.mat');
% Assume the variable is named `ecg`, with 8 rows (leads) and N columns (samples)
Fs = 500;
N = 5000;
t = (0:N-1)/Fs;

% 1. Plot the 8 leads
figure;
for i = 1:8
    subplot(8,1,i);
    plot(t, Y(i,:));
    ylabel(['Lead ', num2str(i)]);
    if i == 1
        title('8-Lead ECG Recording (Leads 1-5: Abdominal, 6-8: Thoracic)');
    end
end
xlabel('Time (s)');

% Observation: Lead 4 (abdominal) appears particularly noisy

% 2. Power Spectral Density (PSD) of the noisy lead (Lead 4)
noisy_lead = Y(4,:);
window = hamming(1024);
nfft = 2048;
[pxx, f] = pwelch(noisy_lead, window, [], nfft, fs);

% Plot PSD
figure;
plot(f, 10*log10(pxx));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of Noisy Abdominal ECG (Lead 4)');
grid on;

% 3. Identify noise frequency intervals from PSD
% Interpretation
fprintf('\nInterpretation of the PSD:\n');
fprintf('- Look for peaks in high-frequency regions (e.g., > 50 Hz) indicating EMG or electronic noise.\n');
fprintf('- Power line interference is often seen around 50 Hz or 60 Hz.\n');
fprintf('- Motion artifacts and baseline wander are seen in low-frequency bands (0.1 - 1 Hz).\n');

% Highlight frequency bands (optional)
hold on;
yL = ylim;
% Highlight power line interference
line([50 50], yL, 'Color', 'r', 'LineStyle', '--');
text(52, yL(2)-10, 'Power Line Noise (50 Hz)', 'Color', 'r');

% Baseline wander band
area(f, 10*log10(pxx), 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3, 'DisplayName', 'Baseline wander region');
xlim([0 100]);  % focus on lower frequencies
legend('PSD', '50 Hz', 'Baseline Wander (~<1 Hz)');


%%

fs = 500;
noisy_lead = Y(4,:);  % 4th abdominal lead

% Welch's method parameters
segmentLength = 512;                % Segment length (you can tune this)
window = hamming(segmentLength);   % Hamming window
overlap = segmentLength / 2;       % 50% overlap
nfft = 2048;                        % Number of FFT points for better resolution

% Estimate PSD using Welch's method
[pxx, f] = pwelch(noisy_lead, window, overlap, nfft, fs);

% Plot PSD
figure;
plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Welch Power Spectral Density - Noisy Abdominal ECG (Lead 4)');
grid on;
xlim([0 100]);  % Focus on ECG + noise band

% Highlight known noise regions
hold on;
yL = ylim;

% Power line interference at 50 Hz
line([50 50], yL, 'Color', 'r', 'LineStyle', '--');
text(52, yL(2)-10, '50 Hz (Power Line)', 'Color', 'r');

% Baseline wander (<1 Hz)
area(f(f < 1), 10*log10(pxx(f < 1)), ...
    'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legend('PSD (Welch)', 'Power Line Noise', 'Baseline Wander Region');

