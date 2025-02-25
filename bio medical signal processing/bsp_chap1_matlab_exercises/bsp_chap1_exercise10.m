% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 10:
%
% Inverse DFT of undersampled Bartlett window's DTFT.
%
%
% HISTORY:
%
% 2025/02/10: - enhanced plots
%
% 2024/02/09: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize = 16;


% generate Bartlett window

N = 64;

w = window(@bartlett, N)';

figure
stem(0:2*N-1, [w, zeros(1, N)])
grid
xlabel('sample value, n')
ylabel('w[n]')


% compute DFT and IDFT

Wdft = fft(w, N);
figure
stem(2*pi*(0:(N-1))/N, abs(Wdft))
grid
hold on
ww = linspace(0, 2*pi, 1000);
W = freqz(w, 1, ww);
plot(ww, abs(W))
set(gca, 'FontSize', fontsize)
xlabel('frequency, \omega (rad/sample)')
ylabel('DTFT, DFT')
title('MSc DSAI - BSP - Spec analysis - exercise 10: 64-point DFT of 64-sample Bartlett window')
legend('W[k]', 'W(e^{j\omega})')

widft = ifft(Wdft);

figure
stem(0:N-1, widft)
grid
set(gca, 'FontSize', fontsize)
title(['IDFT for N = ', num2str(N)]) 
xlabel('sample value, n')
ylabel('w[n]')
title('MSc DSAI - BSP - Spec analysis - exercise 10: 64-point IDFT of 64-sample Bartlett window')


% downsample the DFT by a factor of 2 to get 32-point DFT

Wdft2 = Wdft(1:2:N);

figure
stem(2*pi*(0:(N/2-1))/(N/2), abs(Wdft2))
grid
hold on
plot(ww, abs(W))
set(gca, 'FontSize', fontsize)
xlabel('frequency, \omega (rad/sample)')
ylabel('DTFT, DFT')
title('MSc DSAI - BSP - Spec analysis - exercise 10: 32-point DFT of 64-sample Bartlett window')
legend('W[k]', 'W(e^{j\omega})')

widft2 = ifft(Wdft2);

figure
stem(0:(N/2-1), widft2)
grid
set(gca, 'FontSize', fontsize)
title(['IDFT for N = ', num2str(N/2)])
xlabel('sample value, n')
ylabel('w[n]')
title('MSc DSAI - BSP - Spec analysis - exercise 10: 32-point IDFT of 64-sample Bartlett window')
