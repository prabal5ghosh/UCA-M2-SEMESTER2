% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 11:
%
% Zero padding.
%
%
% HISTORY:
%
% 2025/02/10: - enhanced plots
%
% 2024/02/09: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize = 16;


% DFT of zero-padded sequence

N = 64; % length of the original rectangular window

w = [ones(1, N), zeros(1, N)]; % zero-padded rectangular window

Np = 2*N;

Wdft = fft(w, Np);  % 128-point DFT of the zero-padded 64-point rectangular window

figure
stem(2*pi*(0:(Np-1))/Np, abs(Wdft))
grid
hold on
ww = linspace(0, 2*pi, 1000);
W = freqz(w(1:N), 1, ww);    % DTFT of original 64-point rectangular window
plot(ww, abs(W))
set(gca, 'FontSize', fontsize)
title('MSc DSAI - BSP - Spec analysis - exercise 11: 128-point DFT of 64-sample rectangular window')
xlabel('frequency, \omega (rad/sample)')
ylabel('DTFT, DFT')
legend('W[k]', 'W(e^{j\omega})')


% IDFT of zero-paded sequence

widft = ifft(Wdft);

figure
stem(0:Np-1, widft)
grid
set(gca, 'FontSize', fontsize)
title('MSc DSAI - BSP - Spec analysis - exercise 11: inverse 128-point DFT')
xlabel('sample value, n')
ylabel('w[n]')
