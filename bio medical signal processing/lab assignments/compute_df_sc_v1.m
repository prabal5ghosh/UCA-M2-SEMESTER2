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