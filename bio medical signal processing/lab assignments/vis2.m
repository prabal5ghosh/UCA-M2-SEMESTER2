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

% Optional: Confirm the total number of patients
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
sgtitle('12-Lead ECG for Subject 1');





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


X = Xa_all(7:12,:,:);  % leads V1 to V6
X = permute(X, [3, 2, 1]);  % shape = [75, 15000, 6]

load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indrecur.mat'); indrecur = indrecur - 1;
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\Atrial_prabal_ghosh\indnonrecur.mat'); indnonrecur = indnonrecur - 1;
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

