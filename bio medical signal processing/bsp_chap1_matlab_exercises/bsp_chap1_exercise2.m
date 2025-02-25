% MSc DSAI - Biomedical signal processing
%
% Spectral analysis - exercise 2:
%
% Using MATLAB stem command, plot in the same figure the discrete sinusoids cos(w0n) and 
% cos((w0 + 2pi)n), with w0 = pi/5 rad/sample, for 0 <= n <= 20.
%
%
% HISTORY:
%
% 2024/01/26: - created by Vicente Zarzoso, Université Côte d'Azur, CNRS, I3S.


fontsize = 16;

% generate and plot 1st sinusoid

n = 0:20;           % sample indices to be computed and plotted 
omega0 = pi/5;      % frequency of the first sinusoid
x1 = cos(omega0*n); % first sinusoid

figure
stem(n, x1, 'og')
hold on; grid on;
set(gca, 'FontSize', fontsize)
xlabel('n (samples)')
ylabel('x_1[n], x_2[n]')
title('MSc DSAI - Biomed Sig Proc - Spectral analysis - exercise 2')


% generate and plot 2nd sinusoid

x2 = cos((omega0 + 2*pi)*n);
stem(n, x2, 'xr')

legend('x_1[n]', 'x_2[n]')