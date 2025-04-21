% extract_fecg_mracnc_full.m
% Full MATLAB script to extract fetal ECG (FECG) by cancelling maternal ECG (MECG)
% using the professor's mracnc function (no time dispersion, M = 0).

clear; close all; clc;

% --- 1) Load 8-lead ECG record ----------------------------------------
% Assumes 'pregnancy_ecg.mat' contains variable 'ecg8' of size [N x 8]
load('C:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\ecg_mother.mat');  % => ecg8 (N×8)
ecg8 = Y;  % rename for consistency with rest of script

fs = 500;                % sampling frequency (Hz)
[N, ~] = size(ecg8);     

%% 2) Set MRANC parameters
M  = 0;                  % no time dispersion
mu = 0.01;               % LMS step size

%% 3) Preallocate and run MRANC
FECG = zeros(N,5);
for lead = 1:5
    x = ecg8(:,lead)';      % 1×N primary
    Z = ecg8(:,6:8)';       % 3×N references

    % Call professor’s MRACNC
    s = zeros(1,N);
    W = zeros(3, 1);  % fixed here!
    
    % Build Toeplitz once (but only first row needed since M = 0)
    ZZ = Z;  % size: 3 x N
    
    % Adaptive loop
    for t = 1:N-1
        y_t  = W' * ZZ(:,t);        
        s(t) = x(t) - y_t;             
        W    = W + 2*mu * s(t) * ZZ(:,t);
    end
    
    s(N) = x(N) - W' * ZZ(:,N);

    % Store column
    FECG(:,lead) = s';
end

%% 4) Quick plot for Lead 1 (first 2 sec)
t = (0:N-1)/fs;
samples2 = min(2 * fs, N);

figure;
subplot(2,1,1);
plot(t(1:samples2), ecg8(1:samples2,1));
title('Abdominal Lead 1 (Mixture)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t(1:samples2), FECG(1:samples2,1));
title('Estimated FECG (Lead 1)');
xlabel('Time (s)'); ylabel('Amplitude');

disp(size(ecg8));

%% 5) Save results
save('FECG_mracnc_estimates.mat', 'FECG');

