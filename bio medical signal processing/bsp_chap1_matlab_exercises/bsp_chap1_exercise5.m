% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 5:
%
% Plot the DTFT magnitude of the 64-sample rectangular window.
%
%
% HISTORY:
%
% 2025/29/01: - added figure title and legend
%
% 2024/02/02: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize= 16;

N = 64; % length of the rectangular window

% create frequency axis
w = linspace(-pi, pi, 1000);

% compute DTFT magnitude

W = abs(sin(w*N/2)./sin(w/2));  % use "./" for sample-by-sample division
maxW = max(W);

% plot
figure
plot(w, 20*log10(W))
set(gca, 'FontSize', fontsize)
xlabel('frequency \omega (rad/sample)')
ylabel('|W(e^{j\omega})|')
title('MSc DSAI - Biomed Sig Proc - Spectral analysis - exercise 5')
axis([min(w), max(w), -10, 1.1*max(20*log10(maxW))])
grid
hold on


% using the 'freqz' MATLAB command

Wfreqz = freqz(ones(1, N), 1, w);
plot(w, 20*log10(abs(Wfreqz)), 'r')

legend('direct calculation', '''freqz'' command')

% compute and plot the difference
figure
plot(w, W - abs(Wfreqz))
set(gca, 'FontSize', fontsize)
axis([min(w), max(w), -0.1, 0.1])
grid
ylabel('difference')
title('MSc DSAI - Biomed Sig Proc - Spectral analysis - exercise 5')