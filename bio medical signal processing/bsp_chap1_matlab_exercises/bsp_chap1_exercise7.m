% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 7:
%
% DTFT of common windows.
% 
%
% HISTORY:
%
% 2025/06/02: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize = 16;

N = 25;

w = linspace(-pi, pi, 1000);

figure

wr = window(@rectwin, N);
Wr = freqz(wr, 1, w);
plot(w, 20*log10(abs(Wr)/max(abs(Wr))), 'b', 'LineWidth', 2)

hold on

wt = window(@bartlett, N);
Wt = freqz(wt, 1, w);
plot(w, 20*log10(abs(Wt)/max(abs(Wt))), 'r', 'LineWidth', 2)


wb = window(@blackman, N);
Wb = freqz(wb, 1, w);
plot(w, 20*log10(abs(Wb)/max(abs(Wb))), 'g', 'LineWidth', 2)

axis([-pi, pi, -70, 5])
grid
set(gca, 'FontSize', fontsize')
xlabel('\omega (rad/sample)')
ylabel('|W(e^{j\omega})| (dB)')
title('MSc DSAI - BSP - Spec analysis - exercise 7: DTFT of usual windows')
legend('rectangular', 'Bartlett', 'Blackman')
