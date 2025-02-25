% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 9:
%
% DFT of rectangular window.
% 
%
% HISTORY:
%
% 2025/02/10: - enhanced plots
%
% 2025/01/29: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize = 16;


% compute and plot 64-point DFT of rectangular window

N = 64;

w = ones(1, N);

W = fft(w, N);  % DFT coefficients
omega_k = 2*pi*(0:N-1)/N;

fig_DTFT = figure;
stem(omega_k, abs(W))
hold on

% superimpose DTFT in same frequency interval

omega = linspace(0, 2*pi, 1000);
WW = freqz(w, 1, omega);
plot(omega, abs(WW));
grid
set(gca, 'FontSize', fontsize)
xlabel('\omega (rad/sample)')
ylabel('|W(e^{j\omega})|')
title('MSc DSAI - BSP - Spec analysis - exercise 9: 64-point DFT of 64-sample rectangular window')
legend('W[k], N = 64', 'W(e^{j\omega})')


% compute and plot IDFT

wi = ifft(W);

figure
stem(0:N-1, wi)
grid
set(gca, 'FontSize', fontsize)
xlabel('n')
ylabel('w_i[n]')
title('MSc DSAI - BSP - Spec analysis - exercise 9: 64-point IDFT')

return  % comment out to execute following part


% compute and plot IDFT of oversampled DFT

omega_2N = 2*pi*(0:2*N-1)/(2*N);   % 2N-point sampling of the DTFT
W2 = freqz(w, 1, omega_2N);

figure(fig_DTFT)
stem(omega_2N, abs(W2), 'gx')
legend('W[k], N = 64', 'W(e^{j\omega})', 'W[k], N = 128')

wi2 = ifft(W2);

figure
stem(0:2*N-1, wi2)
grid
set(gca, 'FontSize', fontsize)
xlabel('n')
ylabel('w_i[n]')
title('MSc DSAI - BSP - Spec analysis - exercise 9: 128-point IDFT')
