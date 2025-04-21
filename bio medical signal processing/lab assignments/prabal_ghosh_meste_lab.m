%% STEP 0: Load Data
% Load raw ECG signals from one .mat file
% Load all Rva parts
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Rva1.mat'); 
Xva1 = Rva1;

load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Rva2.mat'); 
Xva2 = Rva2;

load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Rva3.mat'); 
Xva3 = Rva3;

% Merge along the 3rd dimension (patient axis)
Xva_all = cat(3, Xva1, Xva2, Xva3);

fprintf('Total patients loaded: %d\\n', size(Xva_all, 3));

%%
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Ra1.mat'); 
Xa1 = Ra1;

load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Ra2.mat'); 
Xa2 = Ra2;

load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\Ra3.mat'); 
Xa3 = Ra3;

Xa_all = cat(3, Xa1, Xa2, Xa3);
%%


lead = 1;          % Lead I
subject = 1;       % First subject

signal = squeeze(Xva_all(lead, :, subject));  % Get 1D signal for one lead and subject
fs = 256;  % Sampling frequency





%% STEP 1: Visualize 12-lead ECG for this subject
figure;
for i = 1:12
    subplot(6, 2, i);
    plot(squeeze(Rva1(i,:,subject)));
    title(['Lead ' num2str(i)]);
    xlabel('Samples'); ylabel('Amplitude');
end
%sgtitle('12-Lead ECG for Subject 1');





%% R-WAVES 2
% Choose the patient you want to visualize
subject = 1;  % Change this to view other patients (1 to 75)

%% Set up plot
figure('Name', ['12-Lead ECG - Patient ' num2str(subject)], 'NumberTitle', 'off');
for lead = 1:12
    subplot(6, 2, lead);
    plot(Xva_all(lead,:,subject));
    title(['Lead ' num2str(lead)]);
    xlabel('Samples'); ylabel('mV');
    grid on;
end
sgtitle(['12-Lead ECG Visualization - Patient ' num2str(subject)]);


%%
signal = squeeze(Xva_all(1,:,1));  % Lead I, Subject 1
fs = 256;
threshold = max(signal) * 0.5;
min_dist = round(0.2 * fs);

LOCS = [];
PKS = [];
last_peak = -inf;

for i = 2:length(signal)-1
    if signal(i) > threshold && signal(i) > signal(i-1) && signal(i) > signal(i+1)
        if isempty(LOCS) || (i - last_peak) > min_dist
            LOCS = [LOCS i];
            PKS = [PKS signal(i)];
            last_peak = i;
        end
    end
end


figure
plot(signal)
hold on
plot(LOCS, PKS, 'ro')
title('R-peaks Detected - Lead I, Subject 1')
xlabel('Samples'); ylabel('Amplitude (mV)')

%% 3

%LOCS = R-peaks detected earlier
R = [];  % This will hold all QRST segments
fs = 256;
before = 20;
after = min(diff(LOCS)) - 21;  % conservative end point to avoid overlap

% Build the R matrix from Lead I, Subject 1
for k = 2:length(LOCS)-1
    segment = Xva_all(1, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Visualize all segments
figure;
plot(R');  % each row = one QRST
title('Aligned QRST Segments (Lead I)');

% Average QRST waveform
mean_QRST = mean(R);
figure;
plot(mean_QRST, 'r', 'LineWidth', 2);
title('Mean QRST Complex');



%% 4

% Step 1: Do SVD on the QRST segment matrix
[U, S, V] = svd(R);  % R is matrix of all QRST segments

% Step 2: Visualize how much energy is in each singular value
figure;
plot(diag(S), 'o-');
title('Singular Values of QRST Matrix');
xlabel('Component #'); ylabel('Energy');

% Step 3: Create a "ventricular subspace" from first 2 components
M = V(:, 1:2);  % QRS + T-wave basis

% Step 4: Project the first beat onto this basis
a = pinv(M) * R(1,:)';        % Coefficients for reconstruction
reconstructed = M * a;        % QRS + T-wave reconstructed

% Step 5: Subtract QRS from original to get atrial signal
residual = R(1,:)' - reconstructed;

% Step 6: Compare visually
figure;
subplot(3,1,1); plot(R(1,:), 'b'); title('Original QRST Segment');
subplot(3,1,2); plot(reconstructed, 'r'); title('Reconstructed QRS+T');
subplot(3,1,3); plot(residual, 'k'); title('Residual = Atrial Activity');

%%

%% STEP 4: Ventricular Activity Subtraction via SVD
% Assume LOCS = R-peaks detected earlier
R = [];  % This will hold all QRST segments
fs = 256;
before = 20; % Samples before R-peak (~78ms)
after = min(diff(LOCS)) - 21; % Conservative end point to avoid overlap

% Build the R matrix from Lead I, Subject 1
for k = 2:length(LOCS)-1
    segment = Xva_all(1, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R);  % R is matrix of all QRST segments

% Create a "ventricular subspace" from the first 2 components
M = V(:, 1:2);  % QRS + T-wave basis

% Initialize a matrix to store all residuals
residuals = zeros(size(R));  % Same size as R

% Loop through all beats to compute residuals
for k = 1:size(R, 1)
    a = pinv(M) * R(k, :)';        % Project current beat onto QRS basis
    qrs_reconstructed = M * a;    % Reconstruct ventricular activity
    residuals(k, :) = R(k, :)' - qrs_reconstructed;  % Subtract it
end

% Plot all residuals on the same figure
figure;
plot(residuals', 'Color', [0.6 0.6 1]);  % Light blue for all
title('Residuals (Atrial Activity) After SVD-based QRS Subtraction');
xlabel('Samples per Beat'); ylabel('Amplitude (mV)');

% Overlay the mean atrial signal
hold on;
plot(mean(residuals), 'r', 'LineWidth', 2);
legend('All Residuals', 'Mean Atrial Signal');


%%


X = Xa_all(7:12,:,:);  % leads V1 to V6
X = permute(X, [3, 2, 1]);  % shape = [75, 15000, 6]


%%
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indrecur.mat'); indrecur = indrecur - 1;
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indnonrecur.mat'); indnonrecur = indnonrecur - 1;


%%   
% Machine learning for AF recurrence classification (applied to Ra) (student choice for 
%the software, Matlab ?,Tensorflow ?) à Use only V1 to V6 leads ( Ra(7:12,:) ) 

y = zeros(75,1); y(indrecur) = 1;

% Flatten data: 15000×6 → 90000 features per subject
X_flat = reshape(X, 75, []);
X_flat = X_flat';

% Transpose label to match net format
y = y';

% Convert to one-hot (needed by patternnet)
targets = full(ind2vec(y + 1));  % MATLAB likes classes to start at 1

% Create and train the network
net = patternnet([64, 32]);  % Two hidden layers
net.trainParam.epochs = 100;
net.trainParam.showWindow = true;

[net, tr] = train(net, X_flat, targets);

% Evaluate
outputs = net(X_flat);
predictions = vec2ind(outputs) - 1;

accuracy = sum(predictions == y) / length(y);
fprintf('Accuracy: %.2f%%\\n', accuracy * 100);



%%
% 1. Initialize a matrix to store all residuals
residuals = zeros(size(R));  % same size as R

% 2. Loop through all beats
for k = 1:size(R,1)
    a = pinv(M) * R(k,:)';        % project current beat onto QRS basis
    qrs_reconstructed = M * a;    % reconstruct ventricular activity
    residuals(k,:) = R(k,:)' - qrs_reconstructed;  % subtract it
end

% 3. Plot all residuals on the same figure
figure;
plot(residuals', 'Color', [0.6 0.6 1]);  % light blue for all
title('Residuals (Atrial Activity) After SVD-based QRS Subtraction');
xlabel('Samples per Beat'); ylabel('Amplitude (mV)');

% 4. Overlay the mean atrial signal
hold on;
plot(mean(residuals), 'r', 'LineWidth', 2);
legend('All Residuals', 'Mean Atrial Signal');



%% STEP 3: QRST Averaging (aligning beats around R-peaks)
window_before = 20;  % samples before R-peak (~78ms)
window_after = 128;  % samples after R-peak (~0.5s)
minRR = min(diff(LOCS));  % shortest RR interval

segments = [];
for k = 2:length(LOCS)-1
    start_idx = LOCS(k) - window_before;
    end_idx = LOCS(k) + minRR - window_before - 1;
    if start_idx > 0 && end_idx < length(signal)
        segment = signal(start_idx:end_idx);
        segments = [segments; segment];
    end
end

% Plot all QRST segments and their average
figure;
plot(segments', 'b'); hold on;
plot(mean(segments), 'r', 'LineWidth', 2);
title('QRST Aligned Segments and Mean Beat');
xlabel('Samples'); ylabel('Amplitude');

%% STEP 4: Ventricular Activity Subtraction via SVD
R = segments;

[U, S, V] = svd(R);           % SVD on the beat segments
M = V(:,1:2);                 % First 2 singular vectors
a = pinv(M) * R(1,:)';        % Project first beat onto ventricular basis
reconstructed = M * a;       % Reconstruct QRS-T
residual = R(1,:)' - reconstructed;  % Atrial residual

% Display subtraction result
figure;
subplot(3,1,1); plot(R(1,:), 'b'); title('Original QRST Segment');
subplot(3,1,2); plot(reconstructed, 'r'); title('Reconstructed QRS (SVD)');
subplot(3,1,3); plot(residual, 'k'); title('Residual (Atrial Fibrillation)');








%% STEP 4: Ventricular Activity Subtraction via SVD
% Assume LOCS = R-peaks detected earlier
R = [];  % This will hold all QRST segments
fs = 256;
before = 20; % Samples before R-peak (~78ms)
after = min(diff(LOCS)) - 21; % Conservative end point to avoid overlap

% Build the R matrix from Lead I, Subject 1
for k = 2:length(LOCS)-1
    segment = Xva_all(1, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R);  % R is matrix of all QRST segments

% Create a "ventricular subspace" from the first 2 components
M = V(:, 1:2);  % QRS + T-wave basis

% Initialize a matrix to store all residuals
residuals = zeros(size(R));  % Same size as R

% Loop through all beats to compute residuals
for k = 1:size(R, 1)
    a = pinv(M) * R(k, :)';        % Project current beat onto QRS basis
    qrs_reconstructed = M * a;    % Reconstruct ventricular activity
    residuals(k, :) = R(k, :)' - qrs_reconstructed;  % Subtract it
end

% Plot all residuals on the same figure
figure;
plot(residuals', 'Color', [0.6 0.6 1]);  % Light blue for all
title('Residuals (Atrial Activity) After SVD-based QRS Subtraction');
xlabel('Samples per Beat'); ylabel('Amplitude (mV)');

% Overlay the mean atrial signal
hold on;
plot(mean(residuals), 'r', 'LineWidth', 2);
legend('All Residuals', 'Mean Atrial Signal');



%%


%% STEP 4: Ventricular Activity Subtraction via SVD

fs = 256; % Sampling frequency
lead = 1; % First lead (Lead I)

% Extract QRST segments from Xva (raw ECG)
R = []; % This will hold all QRST segments
before = 20; 
after = min(diff(LOCS)) - 21; 

% Build the R matrix from Lead I
for k = 2:length(LOCS)-1
    segment = Xva_all(lead, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R); % R is the matrix of all QRST segments

% Create a "ventricular subspace" from the first 2 components
M = [V(:, 1), V(:, 2)]; % QRS + T-wave basis

% Project the first beat onto the ventricular subspace
a = pinv(M) * R(1, :)'; % Coefficients for reconstruction
reconstructed = M * a; % Reconstruct ventricular activity

% Subtract ventricular activity to get atrial residual
residual = R(1, :)' - reconstructed;

% Plot the result of ventricular activity subtraction
figure;
subplot(3, 1, 1);
plot(R(1, :)', 'b');
title('Original QRST Segment');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 2);
plot(reconstructed, 'r');
title('Reconstructed Ventricular Activity (QRS+T)');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 3);
plot(residual, 'k');
title('Residual Atrial Activity');
xlabel('Samples'); ylabel('Amplitude');

% Extract corresponding atrial activity segment from Xa
Ra = []; % This will hold atrial activity segments
for k = 2:length(LOCS)-1
    segment = Xa_all(lead, LOCS(k)-before : LOCS(k)+after, 1);
    Ra(k-1, :) = segment;
end

% Compare residual atrial activity with Xa
figure;
hold off;
plot(Ra(1, :)', 'g', 'LineWidth', 1.5); % Atrial activity from Xa
hold on;
plot(residual, 'k', 'LineWidth', 1.5); % Residual atrial activity from Xva
title('Comparison of Atrial Activity (Xa vs Residual)');
xlabel('Samples'); ylabel('Amplitude');
legend('Xa (Atrial Activity)', 'Residual (From Xva)');





%%

%% STEP 4: Ventricular Activity Subtraction via SVD
% Assume LOCS = R-peaks detected earlier
fs = 256; % Sampling frequency
lead = 1; % First lead (Lead I)

% Extract QRST segments from Xva (raw ECG)
R = []; % This will hold all QRST segments
before = 20; % Samples before R-peak (~78ms)
after = min(diff(LOCS)) - 21; % Conservative end point to avoid overlap

% Build the R matrix from Lead I
for k = 2:length(LOCS)-1
    segment = Xva_all(lead, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R); % R is the matrix of all QRST segments

% Create a "ventricular subspace" from the first 2 components
M = [V(:, 1), V(:, 2)]; % QRS + T-wave basis

% Loop through all beats to compute residuals
residuals = zeros(size(R)); % Initialize residuals matrix
for k = 1:size(R, 1)
    a = pinv(M) * R(k, :)';        % Project current beat onto QRS basis
    qrs_reconstructed = M * a;    % Reconstruct ventricular activity
    residuals(k, :) = R(k, :)' - qrs_reconstructed;  % Subtract it
end

% Plot the result of ventricular activity subtraction for the first beat
figure;
subplot(3, 1, 1);
plot(R(1, :)', 'b');
title('Original QRST Segment');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 2);
plot(qrs_reconstructed, 'r');
title('Reconstructed Ventricular Activity (QRS+T)');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 3);
plot(residuals(1, :)', 'k');
title('Residual Atrial Activity (After Ventricular Subtraction)');
xlabel('Samples'); ylabel('Amplitude');

% Plot all residuals on the same figure
figure;
plot(residuals', 'Color', [0.6 0.6 1]); % Light blue for all residuals
title('Residuals (Atrial Activity) After SVD-based Ventricular Subtraction');
xlabel('Samples per Beat'); ylabel('Amplitude (mV)');

% Overlay the mean atrial signal
hold on;
plot(mean(residuals), 'r', 'LineWidth', 2);
legend('All Residuals', 'Mean Residual Signal');







%%

% Assume LOCS = R-peaks detected earlier
fs = 256; % Sampling frequency
lead = 1; % First lead (Lead I)

% Extract QRST segments from Xva (raw ECG)
R = []; % This will hold all QRST segments
before = 20; % Samples before R-peak (~78ms)
after = min(diff(LOCS)) - 21; % Conservative end point to avoid overlap

% Build the R matrix from Lead I
for k = 2:length(LOCS)-1
    segment = Xva_all(lead, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R); % R is the matrix of all QRST segments

% Visualize the singular values to understand energy distribution
figure;
plot(diag(S), 'o-');
title('Singular Values of QRST Matrix');
xlabel('Component #'); ylabel('Energy');

% Create a "ventricular subspace" from the first 2 components
M = [V(:, 1), V(:, 2)]; % QRS + T-wave basis

% Project the first beat onto the ventricular subspace
a = pinv(M) * R(1, :)'; % Coefficients for reconstruction
reconstructed = M * a; % Reconstruct ventricular activity

% Subtract ventricular activity to get atrial residual
residual = R(1, :)' - reconstructed;

% Plot the result of ventricular activity subtraction
figure;
subplot(3, 1, 1);
plot(R(1, :)', 'b');
title('Original QRST Segment');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 2);
plot(reconstructed, 'r');
title('Reconstructed Ventricular Activity (QRS+T)');
xlabel('Samples'); ylabel('Amplitude');

subplot(3, 1, 3);
plot(residual, 'k');
title('Residual Atrial Activity');
xlabel('Samples'); ylabel('Amplitude');

% Segmentation of the Xa matrix with the same time location as Xva
Ra = []; % This will hold atrial activity segments
for k = 2:length(LOCS)-1
    segment = Xa_all(lead, LOCS(k)-before : LOCS(k)+after, 1);
    Ra(k-1, :) = segment;
end

% Compare residual atrial activity with Xa
figure;
hold off;
plot(Ra(1, :)', 'g', 'LineWidth', 1.5); % Atrial activity from Xa
hold on;
plot(residual, 'k', 'LineWidth', 1.5); % Residual atrial activity from Xva
title('Comparison of Atrial Activity (Xa vs Residual)');
xlabel('Samples'); ylabel('Amplitude');
legend('Xa (Atrial Activity)', 'Residual (From Xva)');


































%%


% Load the data
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\ecg_mother.mat');
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



%%

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



%%
% High-pass filter design
hp_cutoff = 1; % Cutoff frequency (Hz)
[b_hp, a_hp] = butter(4, hp_cutoff / (fs / 2), 'high'); % 4th-order Butterworth filter

% Apply the high-pass filter
filtered_signal_hp = filtfilt(b_hp, a_hp, noisy_lead);



%%

% Custom Notch Filter Design
notch_freq = 50; % Notch frequency (Hz)
notch_bw = 1; % Bandwidth (Hz)
wo = notch_freq / (fs / 2); % Normalized frequency
bw = notch_bw / (fs / 2); % Normalized bandwidth

% Design the notch filter
[b_notch, a_notch] = butter(2, [wo - bw/2, wo + bw/2], 'stop'); % 2nd-order Butterworth bandstop filter

% Apply the notch filter
filtered_signal_notch = filtfilt(b_notch, a_notch, filtered_signal_hp);





%%

% Low-pass filter design
lp_cutoff = 40; % Cutoff frequency (Hz)
[b_lp, a_lp] = butter(4, lp_cutoff / (fs / 2), 'low'); % 4th-order Butterworth filter

% Apply the low-pass filter
filtered_signal = filtfilt(b_lp, a_lp, filtered_signal_notch);
























%%

% Apply all filters sequentially
filtered_signal_hp = filtfilt(b_hp, a_hp, noisy_lead); % High-pass filter
filtered_signal_notch = filtfilt(b_notch, a_notch, filtered_signal_hp); % Notch filter
filtered_signal = filtfilt(b_lp, a_lp, filtered_signal_notch); % Low-pass filter

% Plot the filtered signal
figure;
plot(t, filtered_signal);
title('Filtered ECG Signal (Lead 4)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');



%%

% Apply all filters sequentially
filtered_signal_hp = filtfilt(b_hp, a_hp, noisy_lead); % High-pass filter
filtered_signal_notch = filtfilt(b_notch, a_notch, filtered_signal_hp); % Notch filter
filtered_signal = filtfilt(b_lp, a_lp, filtered_signal_notch); % Low-pass filter


%%

% Plot original and filtered signals in the time domain
figure;
subplot(2, 1, 1);
plot(t, noisy_lead, 'r');
title('Original Noisy ECG Signal (Lead 4)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

subplot(2, 1, 2);
plot(t, filtered_signal, 'b');
title('Filtered ECG Signal (Lead 4)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;


%%

% Compute PSD of original and filtered signals
[pxx_original, f] = pwelch(noisy_lead, window, overlap, nfft, fs);
[pxx_filtered, f] = pwelch(filtered_signal, window, overlap, nfft, fs);

% Plot PSDs
figure;
plot(f, 10*log10(pxx_original), 'r', 'LineWidth', 1.5); hold on;
plot(f, 10*log10(pxx_filtered), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD of Original and Filtered ECG Signals');
legend('Original Signal', 'Filtered Signal');
grid on;
xlim([0 100]); % Focus on the relevant frequency range



















































%%   

% 1.2 Atrial brillation analysis 


%  1. Create MATLAB functions to compute the DF and SC of the signals given as inputs. As
 %a PSD estimator, use the Welch method with the typical window type and size, overlap
 %factor and number of FFT points found in the literature. Find also in the literature the
 %classical frequency band for SC computation


%% STEP 4: Ventricular Activity Subtraction via SVD
% Assume LOCS = R-peaks detected earlier
R = [];  % This will hold all QRST segments
fs = 256;
before = 20;
after = min(diff(LOCS)) - 21;  % conservative end point to avoid overlap

% Build the R matrix from Lead I, Subject 1
for k = 2:length(LOCS)-1
    segment = Xva_all(1, LOCS(k)-before : LOCS(k)+after, 1);
    R(k-1, :) = segment;
end

% Perform SVD on the QRST segment matrix
[U, S, V] = svd(R);  % R is matrix of all QRST segments

% Create a "ventricular subspace" from the first 2 components
M = V(:, 1:2);  % QRS + T-wave basis

% Initialize a matrix to store all residuals
residuals = zeros(size(R));  % same size as R

% Loop through all beats to compute residuals
for k = 1:size(R, 1)
    a = pinv(M) * R(k, :)';        % Project current beat onto QRS basis
    qrs_reconstructed = M * a;    % Reconstruct ventricular activity
    residuals(k, :) = R(k, :)' - qrs_reconstructed;  % Subtract it
end

% Plot all residuals on the same figure
figure;
plot(residuals', 'Color', [0.6 0.6 1]);  % light blue for all
title('Residuals (Atrial Activity) After SVD-based QRS Subtraction');
xlabel('Samples per Beat'); ylabel('Amplitude (mV)');

% Overlay the mean atrial signal
hold on;
plot(mean(residuals), 'r', 'LineWidth', 2);
legend('All Residuals', 'Mean Atrial Signal');




%%
% 2. Compute the DF and SC in lead V1 of all ECG records in the provided AF dataset.

% filepath: c:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\compute_df_sc_v1.m

% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_results = zeros(num_patients, 1); % Dominant Frequency for each patient
SC_results = zeros(num_patients, 1); % Spectral Concentration for each patient

% Loop through all patients
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal = squeeze(Xa_all(lead, :, patient));
    
    % Compute DF and SC using the provided function
    [DF, SC] = compute_df_sc(signal, fs);
    
    % Store results
    DF_results(patient) = DF;
    SC_results(patient) = SC;
end

% Display summary results
fprintf('Mean Dominant Frequency (DF): %.2f Hz\n', mean(DF_results));
fprintf('Mean Spectral Concentration (SC): %.2f %%\n', mean(SC_results));

% Plot DF and SC distributions
figure;
subplot(2, 1, 1);
histogram(DF_results, 10, 'FaceColor', 'b');
title('Dominant Frequency (DF) Distribution');
xlabel('Frequency (Hz)'); ylabel('Count');

subplot(2, 1, 2);
histogram(SC_results, 10, 'FaceColor', 'r');
title('Spectral Concentration (SC) Distribution');
xlabel('SC (%)'); ylabel('Count');







%%

% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_results = zeros(num_patients, 1); % Dominant Frequency for each patient
SC_results = zeros(num_patients, 1); % Spectral Concentration for each patient

% Loop through all patients
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal = squeeze(Xa_all(lead, :, patient));
    
    % Compute DF and SC using the provided function
    [DF, SC] = compute_df_sc(signal, fs);
    
    % Store results
    DF_results(patient) = DF;
    SC_results(patient) = SC;
end

% Display summary results
fprintf('Mean Dominant Frequency (DF): %.2f Hz\n', mean(DF_results));
fprintf('Mean Spectral Concentration (SC): %.2f %%\n', mean(SC_results));

% Plot DF and SC distributions
figure;
subplot(2, 1, 1);
histogram(DF_results, 10, 'FaceColor', 'b');
title('Dominant Frequency (DF) Distribution');
xlabel('Frequency (Hz)'); ylabel('Count');

subplot(2, 1, 2);
histogram(SC_results, 10, 'FaceColor', 'r');
title('Spectral Concentration (SC) Distribution');
xlabel('SC (%)'); ylabel('Count');



%% 3. With the help of box-and-whiskers plots, compare the distributions of the spectral param
%eters of the recorded data (Xva) versus those of the estimated atrial activity (Xa) over the
 %patient population
% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xva_all, 3); % Total number of patients

% Initialize results
DF_Xva = zeros(num_patients, 1); % Dominant Frequency for Xva
SC_Xva = zeros(num_patients, 1); % Spectral Concentration for Xva
DF_Xa = zeros(num_patients, 1);  % Dominant Frequency for Xa
SC_Xa = zeros(num_patients, 1);  % Spectral Concentration for Xa

% Loop through all patients
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal_Xva = squeeze(Xva_all(lead, :, patient)); % Raw ECG (Xva)
    signal_Xa = squeeze(Xa_all(lead, :, patient));   % Atrial activity (Xa)
    
    % Compute DF and SC for Xva
    [DF_Xva(patient), SC_Xva(patient)] = compute_df_sc(signal_Xva, fs);
    
    % Compute DF and SC for Xa
    [DF_Xa(patient), SC_Xa(patient)] = compute_df_sc(signal_Xa, fs);
end

% Combine data for box-and-whiskers plots
DF_data = [DF_Xva; DF_Xa];
SC_data = [SC_Xva; SC_Xa];
group_labels = categorical([repmat("Xva (Raw ECG)", num_patients, 1); repmat("Xa (Atrial Activity)", num_patients, 1)]);

% Create box-and-whiskers plots
figure;

% Dominant Frequency (DF) comparison
subplot(2, 1, 1);
boxchart(group_labels, DF_data, 'BoxFaceColor', 'b');
title('Comparison of Dominant Frequency (DF)');
ylabel('Frequency (Hz)');

% Spectral Concentration (SC) comparison
subplot(2, 1, 2);
boxchart(group_labels, SC_data, 'BoxFaceColor', 'r');
title('Comparison of Spectral Concentration (SC)');
ylabel('Spectral Concentration (%)');




%% 4. Design a classi er to predict AF recurrence using DF and SC as features. Show its per
%formance using 5-fold cross-validation.

% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xva_all, 3); % Total number of patients

% Initialize results
DF_Xa = zeros(num_patients, 1);  % Dominant Frequency for Xa
SC_Xa = zeros(num_patients, 1);  % Spectral Concentration for Xa

% Extract recurrence labels
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indrecur.mat'); 
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indnonrecur.mat'); 
y = zeros(num_patients, 1); % Initialize labels
y(indrecur) = 1; % Recurrence patients labeled as 1
y(indnonrecur) = 0; % Non-recurrence patients labeled as 0

% Loop through all patients to compute DF and SC for Xa
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal_Xa = squeeze(Xa_all(lead, :, patient)); % Atrial activity (Xa)
    
    % Compute DF and SC for Xa
    [DF_Xa(patient), SC_Xa(patient)] = compute_df_sc(signal_Xa, fs);
end

% Combine DF and SC as features
X_features = [DF_Xa, SC_Xa]; % Feature matrix (DF and SC)

% 5-Fold Cross-Validation
cv = cvpartition(y, 'KFold', 5); % Create 5-fold cross-validation partition
accuracy = zeros(cv.NumTestSets, 1); % Initialize accuracy for each fold

for fold = 1:cv.NumTestSets
    % Split data into training and testing sets
    train_idx = training(cv, fold);
    test_idx = test(cv, fold);
    
    X_train = X_features(train_idx, :);
    y_train = y(train_idx);
    X_test = X_features(test_idx, :);
    y_test = y(test_idx);
    
    % Train Logistic Regression Model
    model = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
    
    % Predict on test set
    y_pred = round(predict(model, X_test)); % Logistic regression outputs probabilities, round to 0 or 1
    
    % Compute accuracy
    accuracy(fold) = sum(y_pred == y_test) / length(y_test);
end

% Display overall performance
mean_accuracy = mean(accuracy) * 100;
fprintf('Mean Accuracy (5-Fold Cross-Validation): %.2f%%\n', mean_accuracy);








%%  5. Repeat parts 2  by averaging, for each patient, the spectral parameters over all
 % precordial leads (V1 to V6).




% Parameters
fs = 256; % Sampling frequency
leads = 7:12; % Precordial leads V1 to V6
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_results = zeros(num_patients, 1); % Averaged Dominant Frequency for each patient
SC_results = zeros(num_patients, 1); % Averaged Spectral Concentration for each patient

% Loop through all patients
for patient = 1:num_patients
    DF_leads = zeros(length(leads), 1); % DF for each lead
    SC_leads = zeros(length(leads), 1); % SC for each lead
    
    % Loop through all precordial leads (V1 to V6)
    for lead_idx = 1:length(leads)
        lead = leads(lead_idx);
        signal = squeeze(Xa_all(lead, :, patient)); % Extract signal for current lead and patient
        
        % Compute DF and SC for the current lead
        [DF_leads(lead_idx), SC_leads(lead_idx)] = compute_df_sc(signal, fs);
    end
    
    % Average DF and SC over all leads (V1 to V6)
    DF_results(patient) = mean(DF_leads);
    SC_results(patient) = mean(SC_leads);
end

% Display summary results
fprintf('Mean Dominant Frequency (DF): %.2f Hz\n', mean(DF_results));
fprintf('Mean Spectral Concentration (SC): %.2f %%\n', mean(SC_results));

% Plot DF and SC distributions
figure;
subplot(2, 1, 1);
histogram(DF_results, 10, 'FaceColor', 'b');
title('Dominant Frequency (DF) Distribution (Averaged Over V1 to V6)');
xlabel('Frequency (Hz)'); ylabel('Count');

subplot(2, 1, 2);
histogram(SC_results, 10, 'FaceColor', 'r');
title('Spectral Concentration (SC) Distribution (Averaged Over V1 to V6)');
xlabel('SC (%)'); ylabel('Count');



%%  5. Repeat parts 4  by averaging, for each patient, the spectral parameters over all
 % precordial leads (V1 to V6).

% Parameters
fs = 256; % Sampling frequency
leads = 7:12; % Precordial leads V1 to V6
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_Xa = zeros(num_patients, 1); % Averaged Dominant Frequency for Xa
SC_Xa = zeros(num_patients, 1); % Averaged Spectral Concentration for Xa

% Extract recurrence labels
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indrecur.mat'); 
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indnonrecur.mat'); 
y = zeros(num_patients, 1); % Initialize labels
y(indrecur) = 1; % Recurrence patients labeled as 1
y(indnonrecur) = 0; % Non-recurrence patients labeled as 0

% Loop through all patients to compute averaged DF and SC for Xa
for patient = 1:num_patients
    DF_leads = zeros(length(leads), 1); % DF for each lead
    SC_leads = zeros(length(leads), 1); % SC for each lead
    
    % Loop through all precordial leads (V1 to V6)
    for lead_idx = 1:length(leads)
        lead = leads(lead_idx);
        signal = squeeze(Xa_all(lead, :, patient)); % Extract signal for current lead and patient
        
        % Compute DF and SC for the current lead
        [DF_leads(lead_idx), SC_leads(lead_idx)] = compute_df_sc(signal, fs);
    end
    
    % Average DF and SC over all leads (V1 to V6)
    DF_Xa(patient) = mean(DF_leads);
    SC_Xa(patient) = mean(SC_leads);
end

% Combine DF and SC as features
X_features = [DF_Xa, SC_Xa]; % Feature matrix (DF and SC)

% Custom 5-Fold Cross-Validation
k = 5; % Number of folds
indices = mod(1:num_patients, k) + 1; % Generate fold indices manually
accuracy = zeros(k, 1); % Initialize accuracy for each fold

for fold = 1:k
    % Split data into training and testing sets
    test_idx = (indices == fold);
    train_idx = ~test_idx;
    
    X_train = X_features(train_idx, :);
    y_train = y(train_idx);
    X_test = X_features(test_idx, :);
    y_test = y(test_idx);
    
    % Train Logistic Regression Model
    model = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
    
    % Predict on test set
    y_pred = round(predict(model, X_test)); % Logistic regression outputs probabilities, round to 0 or 1
    
    % Compute accuracy
    accuracy(fold) = sum(y_pred == y_test) / length(y_test);
end

% Display overall performance
mean_accuracy = mean(accuracy) * 100;
fprintf('Mean Accuracy (5-Fold Cross-Validation): %.2f%%\n', mean_accuracy);



%%
% 6. Conclude on the ability of spectral parameters DF and SC to predict AF recurrence in the
% given dataset.


%  answer- The spectral parameters Dominant Frequency (DF) and Spectral Concentration (SC) show promising potential for predicting AF recurrence in the given dataset. While the achieved classification accuracy is encouraging, further improvements can be made by incorporating additional features, advanced classifiers, and larger datasets. These parameters provide a solid foundation for non-invasive prediction of AF recurrence and could be valuable in clinical applications.
















%%--------------------------------------------------------------------






















%%

function [DF, SC] = compute_df_sc(signal, fs)
    % Function to compute Dominant Frequency (DF) and Spectral Concentration (SC)
    % Inputs:
    %   signal - Input AF signal
    %   fs - Sampling frequency (Hz)
    % Outputs:
    %   DF - Dominant Frequency (Hz)
    %   SC - Spectral Concentration (%)

    % Parameters for Welch's method
    segmentLength = 512;                % Segment length
    window = hamming(segmentLength);    % Hamming window
    overlap = segmentLength / 2;        % 50% overlap
    nfft = 2048;                        % Number of FFT points

    % Compute PSD using Welch's method
    [pxx, f] = pwelch(signal, window, overlap, nfft, fs);

    % Define the AF band for analysis (3–9 Hz)
    af_band = [3, 9]; % AF band in Hz
    af_band_idx = (f >= af_band(1)) & (f <= af_band(2)); % Indices for AF band
    f_af = f(af_band_idx); % Frequencies in AF band
    pxx_af = pxx(af_band_idx); % PSD values in AF band

    % Compute Dominant Frequency (DF)
    [~, max_idx] = max(pxx_af); % Find the peak in the AF band
    DF = f_af(max_idx); % Dominant frequency

    % Compute Spectral Concentration (SC)
    narrow_band_width = 0.5; % Narrow band width around DF (±0.5 Hz)
    narrow_band_idx = (f_af >= (DF - narrow_band_width)) & (f_af <= (DF + narrow_band_width));
    narrow_band_power = sum(pxx_af(narrow_band_idx)); % Power in narrow band
    total_af_power = sum(pxx_af); % Total power in AF band
    SC = (narrow_band_power / total_af_power) * 100; % Spectral concentration as a percentage
end


%%

fs = 256; % Sampling frequency (Hz)
signal = squeeze(Xa_all(7, :, 1)); % Example signal (Lead V1, Patient 1)

% Compute DF and SC
[DF, SC] = compute_df_sc(signal, fs);

% Display results
fprintf('Dominant Frequency (DF): %.2f Hz\n', DF);
fprintf('Spectral Concentration (SC): %.2f %%\n', SC);




%%


% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_results = zeros(num_patients, 1); % Dominant Frequency for each patient
SC_results = zeros(num_patients, 1); % Spectral Concentration for each patient

% Loop through all patients
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal = squeeze(Xa_all(lead, :, patient));
    
    % Compute DF and SC using the provided function
    [DF, SC] = compute_df_sc(signal, fs);
    
    % Store results
    DF_results(patient) = DF;
    SC_results(patient) = SC;
end

% Display summary results
fprintf('Mean Dominant Frequency (DF): %.2f Hz\n', mean(DF_results));
fprintf('Mean Spectral Concentration (SC): %.2f %%\n', mean(SC_results));

% Plot DF and SC distributions
figure;
subplot(2, 1, 1);
histogram(DF_results, 10, 'FaceColor', 'b');
title('Dominant Frequency (DF) Distribution');
xlabel('Frequency (Hz)'); ylabel('Count');

subplot(2, 1, 2);
histogram(SC_results, 10, 'FaceColor', 'r');
title('Spectral Concentration (SC) Distribution');
xlabel('SC (%)'); ylabel('Count');



%%


% Parameters
fs = 256; % Sampling frequency
lead = 7; % Lead V1 (7th lead in the dataset)
num_patients = size(Xva_all, 3); % Total number of patients

% Initialize results
DF_Xva = zeros(num_patients, 1); % Dominant Frequency for Xva
SC_Xva = zeros(num_patients, 1); % Spectral Concentration for Xva
DF_Xa = zeros(num_patients, 1);  % Dominant Frequency for Xa
SC_Xa = zeros(num_patients, 1);  % Spectral Concentration for Xa

% Loop through all patients
for patient = 1:num_patients
    % Extract the signal for lead V1 of the current patient
    signal_Xva = squeeze(Xva_all(lead, :, patient)); % Raw ECG (Xva)
    signal_Xa = squeeze(Xa_all(lead, :, patient));   % Atrial activity (Xa)
    
    % Compute DF and SC for Xva
    [DF_Xva(patient), SC_Xva(patient)] = compute_df_sc(signal_Xva, fs);
    
    % Compute DF and SC for Xa
    [DF_Xa(patient), SC_Xa(patient)] = compute_df_sc(signal_Xa, fs);
end

% Combine data for box-and-whiskers plots
DF_data = [DF_Xva; DF_Xa];
SC_data = [SC_Xva; SC_Xa];
group_labels = categorical([repmat("Xva (Raw ECG)", num_patients, 1); repmat("Xa (Atrial Activity)", num_patients, 1)]);

% Create box-and-whiskers plots
figure;

% Dominant Frequency (DF) comparison
subplot(2, 1, 1);
boxchart(group_labels, DF_data, 'BoxFaceColor', 'b');
title('Comparison of Dominant Frequency (DF)');
ylabel('Frequency (Hz)');

% Spectral Concentration (SC) comparison
subplot(2, 1, 2);
boxchart(group_labels, SC_data, 'BoxFaceColor', 'r');
title('Comparison of Spectral Concentration (SC)');
ylabel('Spectral Concentration (%)');



%%

% Parameters
fs = 256; % Sampling frequency
leads = 7:12; % Precordial leads V1 to V6
num_patients = size(Xa_all, 3); % Total number of patients

% Initialize results
DF_results = zeros(num_patients, 1); % Averaged Dominant Frequency for each patient
SC_results = zeros(num_patients, 1); % Averaged Spectral Concentration for each patient

% Loop through all patients
for patient = 1:num_patients
    DF_leads = zeros(length(leads), 1); % DF for each lead
    SC_leads = zeros(length(leads), 1); % SC for each lead
    
    % Loop through all precordial leads (V1 to V6)
    for lead_idx = 1:length(leads)
        lead = leads(lead_idx);
        signal = squeeze(Xa_all(lead, :, patient)); % Extract signal for current lead and patient
        
        % Compute DF and SC for the current lead
        [DF_leads(lead_idx), SC_leads(lead_idx)] = compute_df_sc(signal, fs);
    end
    
    % Average DF and SC over all leads (V1 to V6)
    DF_results(patient) = mean(DF_leads);
    SC_results(patient) = mean(SC_leads);
end

% Display summary results
fprintf('Mean Dominant Frequency (DF): %.2f Hz\n', mean(DF_results));
fprintf('Mean Spectral Concentration (SC): %.2f %%\n', mean(SC_results));

% Plot DF and SC distributions
figure;
subplot(2, 1, 1);
histogram(DF_results, 10, 'FaceColor', 'b');
title('Dominant Frequency (DF) Distribution (Averaged Over V1 to V6)');
xlabel('Frequency (Hz)'); ylabel('Count');

subplot(2, 1, 2);
histogram(SC_results, 10, 'FaceColor', 'r');
title('Spectral Concentration (SC) Distribution (Averaged Over V1 to V6)');
xlabel('SC (%)'); ylabel('Count');




%%


% Extract recurrence labels
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indrecur.mat'); 
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indnonrecur.mat'); 
y = zeros(num_patients, 1); % Initialize labels
y(indrecur) = 1; % Recurrence patients labeled as 1
y(indnonrecur) = 0; % Non-recurrence patients labeled as 0

% Combine DF and SC as features
X_features = [DF_results, SC_results]; % Feature matrix (averaged DF and SC)

% Custom 5-Fold Cross-Validation
k = 5; % Number of folds
indices = mod(1:num_patients, k) + 1; % Generate fold indices manually
accuracy = zeros(k, 1); % Initialize accuracy for each fold

for fold = 1:k
    % Split data into training and testing sets
    test_idx = (indices == fold);
    train_idx = ~test_idx;
    
    X_train = X_features(train_idx, :);
    y_train = y(train_idx);
    X_test = X_features(test_idx, :);
    y_test = y(test_idx);
    
    % Train Logistic Regression Model
    model = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
    
    % Predict on test set
    y_pred = round(predict(model, X_test)); % Logistic regression outputs probabilities, round to 0 or 1
    
    % Compute accuracy
    accuracy(fold) = sum(y_pred == y_test) / length(y_test);
end

% Display overall performance
mean_accuracy = mean(accuracy) * 100;
fprintf('Mean Accuracy (5-Fold Cross-Validation): %.2f%%\n', mean_accuracy);



%%









































































%%

% Generate synthetic signals based on Exercise 6
fs = 500; % Sampling frequency
N = 1000; % Number of samples
t = (0:N-1)/fs; % Time vector

% Primary signal: Desired signal + Interference
f_desired = 5; % Desired signal frequency (5 Hz)
f_interference = 50; % Interference frequency (50 Hz)

% Create desired signal (clean signal)
d = sin(2*pi*f_desired*t);

% Create reference signals (pure interference)
z1 = sin(2*pi*f_interference*t);
z2 = cos(2*pi*f_interference*t);
Z = [z1; z2]; % Reference signals

% Parameter combinations to test
mu_values = [0.0005, 0.005, 0.01, 0.01, 0.01];
M_values = [0, 0, 5, 50, 550];

% Initialize results table
results = zeros(length(mu_values), 3); % [SIR_out, SIR_improvement, convergence_time]

% Test all parameter combinations
for i = 1:length(mu_values)
    mu = mu_values(i);
    M = M_values(i);
    
    % Create interference with controlled SIR
    SIR_in = 0; % Input SIR in dB
    A_interference = 10^(-SIR_in/20); % Interference amplitude
    interference = A_interference * sin(2*pi*f_interference*t);
    
    % Primary input: desired signal + interference
    x = d + interference;
    
    % Apply MRANC algorithm
    [s, W] = mranc(x, Z, M, mu);
    
    % Compute performance metrics
    % 1. Wait for convergence (find steady-state)
    conv_start = round(0.7*N); % Use last 30% of samples for steady-state analysis
    
    % 2. Calculate output SIR
    residual_noise = s(conv_start:end) - d(conv_start:end);
    SIR_out = 10*log10(var(d(conv_start:end))/var(residual_noise));
    
    % 3. Calculate SIR improvement
    original_noise = x(conv_start:end) - d(conv_start:end);
    SIR_in_actual = 10*log10(var(d(conv_start:end))/var(original_noise));
    SIR_improvement = SIR_out - SIR_in_actual;
    
    % 4. Estimate convergence time
    error_curve = abs(s - d);
    error_smooth = movmean(error_curve, 50); % Smooth error curve for stable detection
    threshold = 1.5 * mean(error_smooth(conv_start:end)); % Threshold at 1.5x steady-state error
    convergence_time = find(error_smooth < threshold, 1);
    % 4. Estimate convergence time
    error_curve = abs(s - d);
    error_smooth = movmean(error_curve, 50); % Smooth error curve for stable detection
    threshold = 1.5 * mean(error_smooth(conv_start:end)); % Threshold at 1.5x steady-state error
    convergence_indices = find(error_smooth < threshold, 1);
    
    % Handle case where no convergence point is found
    if isempty(convergence_indices)
        convergence_time = NaN; % or could use length(s) as default
    else
        convergence_time = convergence_indices;
    end
    % Store results
    results(i, :) = [SIR_out, SIR_improvement, convergence_time];
    
    % Plot relevant signals for this parameter combination
    figure;
    
    % Time-domain signals
    subplot(3,1,1);
    plot(t, d, 'b', t, x, 'r', t, s, 'g', 'LineWidth', 1.5);
    title(sprintf('MRANC with μ = %.4f, M = %d', mu, M));
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Desired Signal', 'Noisy Input', 'MRANC Output');
    grid on;
    
    % Error signal
    subplot(3,1,2);
    plot(t, error_curve, 'r', t, error_smooth, 'b', 'LineWidth', 1.5);
    if ~isnan(convergence_time)
        hold on;
        plot([t(convergence_time), t(convergence_time)], ylim, 'k--');
        text(t(convergence_time)+0.1, 0, sprintf('Convergence: %d samples', convergence_time));
    end
    title('Error Signal |e(t)| = |s(t) - d(t)|');
    xlabel('Time (s)');
    ylabel('Error Magnitude');
    grid on;
    
    % Frequency spectrum
    subplot(3,1,3);
    [Px, f] = pwelch(x, hamming(256), 128, 1024, fs);
    [Ps, f] = pwelch(s, hamming(256), 128, 1024, fs);
    plot(f, 10*log10(Px), 'r', f, 10*log10(Ps), 'g', 'LineWidth', 1.5);
    title('Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Noisy Input', 'MRANC Output');
    grid on;
    xlim([0, 100]);
end

% Display results in a table
fprintf('-----------------------------------------------------------------------\n');
fprintf('| μ      | M     | Output SIR (dB) | SIR Improvement (dB) | Convergence Time |\n');
fprintf('-----------------------------------------------------------------------\n');
for i = 1:length(mu_values)
    if isnan(results(i, 3))
        fprintf('| %.4f | %5d | %13.2f | %19.2f | %15s |\n', ...
            mu_values(i), M_values(i), results(i, 1), results(i, 2), 'Did not converge');
    else
        fprintf('| %.4f | %5d | %13.2f | %19.2f | %15d |\n', ...
            mu_values(i), M_values(i), results(i, 1), results(i, 2), results(i, 3));
    end
end
fprintf('-----------------------------------------------------------------------\n');



%%

N = 1000;
n = 0:N-1;
s = sin(2*pi*0.05*n);                % desired signal
noise = randn(1, N);                 % primary noise
z1 = filter([1 0.5], 1, noise);      % reference noise 1
z2 = filter([1 -0.3], 1, noise);     % reference noise 2
Z = [z1; z2];                        % multi-reference matrix
x = s + 0.8*noise;                   % primary signal with noise



mu = 0.005; M = 0;
[s_out, W] = mracnc(x, Z, M, mu);


% After steady-state: use last 500 samples, for example
clean = s(end-499:end);
noisy = x(end-499:end);
residual = s_out(end-499:end);

SIR_before = 10*log10(sum(clean.^2) / sum((noisy - clean).^2));
SIR_after = 10*log10(sum(clean.^2) / sum((residual).^2));
SIR_improvement = SIR_after - SIR_before;



plot(n, x, 'r'); hold on;
plot(n, s_out, 'b');
legend('Noisy signal', 'Filtered output');
xlabel('Samples'); ylabel('Amplitude');
title(['MRANC Output (μ = ', num2str(mu), ', M = ', num2str(M), ')']);




%%


% Generate synthetic signals
fs = 500; % Sampling frequency
N = 1000; % Number of samples
t = (0:N-1)/fs; % Time vector

% Desired signal: unit-amplitude sinusoid with frequency 5 Hz
f_desired = 5;
d = sin(2*pi*f_desired*t);

% Reference signals: sine and cosine at 50 Hz (power line interference)
f_interference = 50;
z1 = sin(2*pi*f_interference*t);
z2 = cos(2*pi*f_interference*t);
Z = [z1; z2]; % Reference signals

% Parameter combinations to test
mu_values = [0.0005, 0.005, 0.01, 0.01, 0.01];
M_values = [0, 0, 5, 50, 550];

% Initialize results table
results = zeros(length(mu_values), 3); % [SIR_out, SIR_improvement, convergence_time]

% Test all parameter combinations
for i = 1:length(mu_values)
    mu = mu_values(i);
    M = M_values(i);
    
    % Create interference with controlled SIR
    SIR_in = 0; % Input SIR in dB (0dB means signal and noise have equal power)
    A_interference = 10^(-SIR_in/20); % Interference amplitude
    interference = A_interference * sin(2*pi*f_interference*t);
    
    % Primary input: desired signal + interference
    x = d + interference;
    
    % Apply MRANC algorithm
    [s, W] = mranc(x, Z, M, mu);
    
    % Compute performance metrics
    % 1. Wait for convergence (find steady-state)
    conv_start = round(0.7*N); % Use last 30% of samples for steady-state analysis
    
    % 2. Calculate output SIR
    residual_noise = s(conv_start:end) - d(conv_start:end);
    SIR_out = 10*log10(var(d(conv_start:end))/var(residual_noise));
    
    % 3. Calculate SIR improvement
    original_noise = x(conv_start:end) - d(conv_start:end);
    SIR_in_actual = 10*log10(var(d(conv_start:end))/var(original_noise));
    SIR_improvement = SIR_out - SIR_in_actual;
    
    % 4. Estimate convergence time
    error_curve = abs(s - d);
    error_smooth = movmean(error_curve, 50); % Smooth error curve for stable detection
    threshold = 1.5 * mean(error_smooth(conv_start:end)); % Threshold at 1.5x steady-state error
    convergence_indices = find(error_smooth < threshold, 1);
    
    % Handle case where no convergence point is found
    if isempty(convergence_indices)
        convergence_time = NaN; % or could use length(s) as default
    else
        convergence_time = convergence_indices;
    end
    
    % Store results
    results(i, :) = [SIR_out, SIR_improvement, convergence_time];
    
    % Plot relevant signals for this parameter combination
    figure;
    
    % Time-domain signals
    subplot(3,1,1);
    plot(t, d, 'b', t, x, 'r', t, s, 'g', 'LineWidth', 1.5);
    title(sprintf('MRANC with μ = %.4f, M = %d', mu, M));
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Desired Signal', 'Noisy Input', 'MRANC Output');
    grid on;
    
    % Error signal
    subplot(3,1,2);
    plot(t, error_curve, 'r', t, error_smooth, 'b', 'LineWidth', 1.5);
    if ~isnan(convergence_time)
        hold on;
        plot([t(convergence_time), t(convergence_time)], ylim, 'k--');
        text(t(convergence_time)+0.1, 0, sprintf('Convergence: %d samples', convergence_time));
    end
    title('Error Signal |e(t)| = |s(t) - d(t)|');
    xlabel('Time (s)');
    ylabel('Error Magnitude');
    grid on;
    
   % Frequency spectrum
subplot(3,1,3);
% Check if s has valid values for PSD calculation
if length(s) >= 256 && ~isempty(s) && all(isfinite(s))
    [Px, f] = pwelch(x, hamming(min(256, length(x))), 128, 1024, fs);
    [Ps, f] = pwelch(s, hamming(min(256, length(s))), 128, 1024, fs);
    plot(f, 10*log10(Px), 'r', f, 10*log10(Ps), 'g', 'LineWidth', 1.5);
    title('Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Noisy Input', 'MRANC Output');
    grid on;
    xlim([0, 100]);
else
    % Just show a text message if PSD can't be calculated
    text(0.5, 0.5, 'Cannot calculate PSD - insufficient data', ...
         'HorizontalAlignment', 'center');
    axis off;
end
end

% Display results in a table
fprintf('-----------------------------------------------------------------------\n');
fprintf('| μ      | M     | Output SIR (dB) | SIR Improvement (dB) | Convergence Time |\n');
fprintf('-----------------------------------------------------------------------\n');
for i = 1:length(mu_values)
    if isnan(results(i, 3))
        fprintf('| %.4f | %5d | %13.2f | %19.2f | %15s |\n', ...
            mu_values(i), M_values(i), results(i, 1), results(i, 2), 'Did not converge');
    else
        fprintf('| %.4f | %5d | %13.2f | %19.2f | %15d |\n', ...
            mu_values(i), M_values(i), results(i, 1), results(i, 2), results(i, 3));
    end
end
fprintf('-----------------------------------------------------------------------\n');

































%%


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


% filepath: c:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\prabal_ghosh_meste_lab.m

% Parameters
fs = 256; % Sampling frequency
leads = 7:12; % Precordial leads V1 to V6
patient = 1; % First patient

% Extract the signals for the first patient
signals = squeeze(Xva_all(leads, :, patient)); % Shape: [6, time]

% Time vector
N = size(signals, 2); % Number of samples
t = (0:N-1) / fs; % Time in seconds

% Plot the 6-lead record
figure;
for i = 1:6
    subplot(6, 1, i);
    plot(t, signals(i, :));
    title(['Lead V', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
end
sgtitle('6-Lead ECG Recording (Precordial Leads V1 to V6) - Patient 1');

% Highlight that atrial activity is more clearly present in lead V1
figure;
plot(t, signals(1, :), 'b');
title('Lead V1 - Atrial Activity');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Zoom into a specific time range (e.g., 0 to 2 seconds)
zoom_start = 0; % Start time in seconds
zoom_end = 2;   % End time in seconds
xlim([zoom_start, zoom_end]);







%%


% filepath: c:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\prabal_ghosh_meste_lab.m

% Parameters
fs = 256; % Sampling frequency
primary_lead = 7; % Lead V1 (primary input)
reference_leads = 8:12; % Leads V2 to V6 (reference signals)
patient = 1; % First patient
mu = 0.01; % Step size for adaptive filter
M = 0; % No time dispersion

% Extract signals for the first patient
primary_signal = squeeze(Xva_all(primary_lead, :, patient)); % Lead V1
reference_signals = squeeze(Xva_all(reference_leads, :, patient))'; % Leads V2-V6 (transpose for proper shape)

% Adaptive Noise Cancellation (ANC) function
function [output_signal, weights] = adaptive_cancellation(primary, references, mu)
    [N, R] = size(references); % N: number of samples, R: number of reference signals
    weights = zeros(R, 1); % Initialize adaptive weights
    output_signal = zeros(N, 1); % Initialize output signal
    for n = 1:N
        ref_vector = references(n, :)'; % Current reference vector (R x 1)
        estimated_qrst = weights' * ref_vector; % Estimated QRST complex
        output_signal(n) = primary(n) - estimated_qrst; % Error signal (residual atrial activity)
        weights = weights + 2 * mu * output_signal(n) * ref_vector; % Update weights
    end
end

% Perform adaptive cancellation
[residual_signal, final_weights] = adaptive_cancellation(primary_signal, reference_signals, mu);

% Time vector
N = length(primary_signal);
t = (0:N-1) / fs;

% Plot results
figure;

% Plot primary signal (Lead V1)
subplot(3, 1, 1);
plot(t, primary_signal, 'r');
title('Primary Signal (Lead V1 - Mixed QRST + Atrial Activity)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot one reference signal (e.g., Lead V2)
subplot(3, 1, 2);
plot(t, reference_signals(:, 1), 'k');
title('Reference Signal (Lead V2 - QRST Complex)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot residual signal (Atrial activity after QRST cancellation)
subplot(3, 1, 3);
plot(t, residual_signal, 'b');
title('Residual Signal (Atrial Activity)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;









%%


% ...existing code...
fs = 500;
N = 1000;
t = (0:N-1)/fs;

% Reduce desired signal power to 0.01 (amplitude = sqrt(0.01) = 0.1)
desired_amplitude = 0.1;
f_desired = 5;  % 5 Hz
d = desired_amplitude * sin(2*pi*f_desired*t);

% ...generate interference signals, mix them, then apply PCA/ICA...

% Perform PCA using bss_pca.m
[W_pca, S_pca] = bss_pca(mixed_signals);

% Perform ICA (e.g., using RobustICA)
% robustica_path = 'path_to_robustica';
% addpath(robustica_path);
% S_ica = robustica(mixed_signals);

% ...existing code to evaluate separation performance...