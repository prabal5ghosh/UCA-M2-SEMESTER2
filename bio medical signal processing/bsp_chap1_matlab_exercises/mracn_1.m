% filepath: /path/to/synthetic_signal.m

% Parameters
N = 1000; % Length of the signal
R = 2; % Number of reference signals
M = 0; % Filter order
mu = 0.01; % Adaptation coefficient
fs = 500; % Sampling frequency
t = (0:N-1)/fs; % Time vector

% Desired signal: unit-amplitude sinusoid with frequency 5 Hz
f_desired = 5;
phase_desired = rand * 2 * pi;
d = sin(2 * pi * f_desired * t + phase_desired);

% Primary signal: desired signal contaminated by a sinusoid of 50 Hz
f_interference = 50;
phase_interference = rand * 2 * pi;
interference = sin(2 * pi * f_interference * t + phase_interference);
SIR = 10; % Signal-to-interference ratio in dB
A_interference = 10^(-SIR/20);
x = d + A_interference * interference;

% Reference signals: unit-amplitude sinus and cosinus with the same frequency as the mains interference
z1 = sin(2 * pi * f_interference * t);
z2 = cos(2 * pi * f_interference * t);
Z = [z1; z2];

% Apply MRANC algorithm
[y, W] = mranc(x, Z, M, mu);

% Plot results
figure;
subplot(4, 1, 1);
plot(t, d);
title('Desired Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4, 1, 2);
plot(t, x);
title('Primary Input');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4, 1, 3);
plot(t, y);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4, 1, 4);
plot(t, d - y);
title('Error Signal');
xlabel('Time (s)');
ylabel('Amplitude');