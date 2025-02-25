% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 6:
%
% Plot the DTFT of a sum of sinusoids with different frequencies and
% amplitudes, for finite observation length.
%
%
% HISTORY:
%
% 2025/29/01: - added figure title
%
% 2024/02/01: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


w = linspace(0, pi, 500);   % plot only positive frequency range

w1 = 0.4*pi;
w2 = 0.45*pi;

A = [1,   1,  1, 0.2];   % 2nd sinusoid amplitude values
N = [64, 32, 40,  40]; % sample size values

for k = 1:length(A)
    x = sumsin(w1, w2, A(k), N(k));
    plot_dtft(x, w);
    title(['MSc DSAI - BSP - Spec analysis - exercise 6: A = ', num2str(A(k)), ', N = ', num2str(N(k))])
    hold on;
    plot(w1*[1, 1], get(gca, 'YLim'), '--k')   % mark sinusoids' frequency locations
    plot(w2*[1, 1], get(gca, 'YLim'), '--k')
    pause
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = sumsin(w1, w2, A, N)
% Generates sum of sinusoids with given frequencies and sample size

    n = 0:N-1;
    x = cos(w1*n) + A*cos(w2*n);
end


function plot_dtft(x, w)
% Computes and plots discrete-time Fourier transform of input sequence

    H = freqz(x, 1, w);  % compute DTFT
    maxabsH = max(abs(H)); maxabsH_dB = 20*log10(maxabsH);
    figure;
    plot(w, 20*log10(abs(H)), 'LineWidth', 2)
    set(gca, 'FontSize', 16); 
    axis([min(w), max(w), -10, 1.1*maxabsH_dB])
    grid
%     set(gca, 'XTick', -pi:pi/2:pi)
%     set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel('\omega (rad/sample)')
    ylabel('|X(e^{j\omega})| (dB)')
end

