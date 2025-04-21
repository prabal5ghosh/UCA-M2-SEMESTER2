load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\ecg_mother.mat');
Fs = 500;
N = 5000;
t = (0:N-1)/Fs;


%%
figure;
for lead = 1:8
    subplot(8,1,lead);
    plot(t, Y(lead, :));
    ylabel(['Lead ', num2str(lead)]);
end
xlabel('Time (s)');




%%
% Calculate PSD of lead4
lead4 = Y(4, :);
[pxx, f] = periodogram(lead4, rectwin(N), N, Fs);

% Calculate the mean of the PSD
mean_psd = mean(10*log10(pxx));

% Plot the PSD
figure;
plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.5); % PSD in blue
hold on;

% Plot the mean of the PSD
yline(mean_psd, 'r--', 'LineWidth', 1.5); % Mean PSD as a red dashed line
hold off;

% Set plot limits and labels
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Power (dB/Hz)');
title('PSD of Lead 4 with Mean Power');
legend('PSD', 'Mean Power');



%%

% Load the data
load('ecg_mother.mat'); % Ensure the variable containing the signal is loaded
Fs = 500; % Sampling frequency in Hz

% Check the size of ecg_mother and adjust N
N = 5000;
signal = ecg_mother(1:N); % Extract the first N samples

% Proceed with PSD estimation
[pxx, f] = periodogram(signal, rectwin(length(signal)), length(signal), Fs);

% Plot the PSD
figure;
plot(f, 10*log10(pxx), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB/Hz)');
title('Estimated PSD of ECG Signal');
grid on;