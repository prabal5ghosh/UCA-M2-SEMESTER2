% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 16:
%
% Spectral analysis of noisy sum of sinusoids via the periodogram.
%
%
% HISTORY:
%
% 2025/02/10: - enhanced plots
%
% 2024/02/07: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.

fontsize = 16;


% generate noisless sum of sinusoids

disp('>>> Generating noiseless sum of sinusoids...')

F1 = 3e2;
F2 = 4e2;   % sinusoid frequencies
T = 2;      % sampling duration in seconds
Fs = 2e3;   % sampling rate
N = T*Fs;   % total number of samples
n = 0:N-1;  % sample values
Nplot = 100;    % number of samples to plot

x = sin(2*pi*F1*n/Fs + 2*pi*rand) + sin(2*pi*F2*n/Fs + 2*pi*rand);

figure
stem(0:Nplot-1, x(1:Nplot))
set(gca, 'FontSize', fontsize)
grid
title('MSc DSAI - BSP - Spec analysis - exercise 16: noiseless sum of sinusoids')
xlabel('sample index, n')
hold on

sound(x, Fs)    % listen to the signal 
pause


% generate noisy observation

disp('>>> Generating noisy observation...')

SNR_dB = 10;        % signal-to-noise ratio in dB
Px = 1/N*sum(x.^2); % power of the noiseless signal
b = sqrt(Px*10^(-SNR_dB/10))*randn(1, N);
y = x + b;

disp(['Generated SNR = ', num2str(10*log10(Px/mean(b.^2))), ' dB'])
stem(0:Nplot-1, y(1:Nplot), 'rx')
title('MSc DSAI - BSP - Spec analysis - exercise 16: noisy sum of sinusoids')
legend('noiseless signal, x[n]', 'noisy signal, y[n]')

sound(y, Fs)
pause


% compute periodogram

disp('>>> Computing the periodogram...')

Nfft = 2^(ceil(log2(N)));   % FFT size

Y = abs(fft(y, Nfft)).^2/N;
w = linspace(0, Fs/2, Nfft/2);

if ~exist('fig_per')
    fig_per = figure;
else
    figure(fig_per);
end

plot(w, 10*log10(Y(1:Nfft/2)))
hold on
axis([0, Fs/2, -30, 40])
grid
ylabel('S_y(\omega)')
xlabel('frequency (Hz)')


return   %%% comment out to do filtering


%%% filtering

disp('>>> Applying lowpass filter...')

load('lowpass_chevyshev_type_II.mat')
z = filter(Hd, y);

sound(z, Fs)

Z = abs(fft(z, Nfft)).^2/N;

if ~exist('fig_per_filt')
    fig_per_filt = figure;
else
    figure(fig_per_filt);
end

plot(w, 10*log10(Z(1:Nfft/2)))
hold on
axis([0, Fs/2, -30, 40])
grid
set(gca, 'FontSize', fontsize)
title('MSc DSAI - BSP - Spec analysis - exercise 16: periodogram')
ylabel('S_z(\omega)')
xlabel('frequency (Hz)')

