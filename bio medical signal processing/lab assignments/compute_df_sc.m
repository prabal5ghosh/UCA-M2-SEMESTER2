% filepath: c:\Users\praba\Documents\GitHub\UCA-M2-SEMESTER2\bio medical signal processing\lab assignments\compute_df_sc.m

function [DF, SC] = compute_df_sc(signal, fs)
    % Compute Dominant Frequency (DF) and Spectral Concentration (SC)
    % Inputs:
    %   signal - Input ECG signal (1D array)
    %   fs - Sampling frequency (Hz)
    % Outputs:
    %   DF - Dominant Frequency (Hz)
    %   SC - Spectral Concentration (%)

    % Parameters for Welch's method
    window = hamming(256); % Typical window size
    overlap = 128;         % 50% overlap
    nfft = 1024;           % Number of FFT points

    % Compute PSD using Welch's method
    [Pxx, f] = pwelch(signal, window, overlap, nfft, fs);

    % Define AF frequency band [3, 9] Hz
    af_band = (f >= 3 & f <= 9);
    f_af = f(af_band);
    Pxx_af = Pxx(af_band);

    % Compute Dominant Frequency (DF)
    [~, peak_idx] = max(Pxx_af);
    DF = f_af(peak_idx);

    % Compute Spectral Concentration (SC)
    narrow_band = (f >= (DF - 0.5) & f <= (DF + 0.5)); % Narrow band around DF
    Pxx_narrow = Pxx(narrow_band);
    SC = (sum(Pxx_narrow) / sum(Pxx_af)) * 100; % Percentage of power in narrow band
end