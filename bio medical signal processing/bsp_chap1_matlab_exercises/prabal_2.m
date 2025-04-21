% complete_lab.m
% Exercises 5 & 6 and Section 2.1 Synthetic Signals (MRANC)

clc; clear; close all;

%% Exercise 5: DTFT of 64-point rectangular window
% Parameters
Nw = 64;
M = Nw - 1;
w = ones(1, Nw);

% Frequency grid
omega = linspace(-pi, pi, 1000);

% Compute DTFT via definition
W_dtft = zeros(size(omega));
for k = 1:length(omega)
    W_dtft(k) = sum(w .* exp(-1j * omega(k) * (0:M)));
end

% Theoretical expression: W(ω) = e^{-jωM/2} * sin(ωNw/2)/sin(ω/2)
W_theo = exp(-1j * omega * M/2) .* (sin(omega * Nw/2) ./ sin(omega/2));
W_theo(omega==0) = Nw;  % handle ω=0

% Using freqz
[H_freqz, omega_f] = freqz(w, 1, 1000, 'whole');
omega_f = omega_f - 2*pi*(omega_f>pi);

% Plot magnitudes
figure;
plot(omega, abs(W_dtft), 'b', omega, abs(W_theo), 'r--', omega_f, abs(H_freqz), 'g:','LineWidth',1.2);
xlabel('Frequency (rad/sample)');
ylabel('|W(e^{jω})|');
legend('DTFT sum', 'Theoretical', 'freqz');
title('DTFT Magnitude of 64-point Rectangular Window');
grid on;

%% Exercise 6: DTFT of length-N sequence x[n] = cos(0.4π n) + A cos(0.45π n)
cases = {
    struct('A',1,'N',64),
    struct('A',1,'N',32),
    struct('A',1,'N',40),
    struct('A',0.2,'N',40)
};
omega = linspace(-pi, pi, 1000);

for idx = 1:length(cases)
    A = cases{idx}.A;
    N = cases{idx}.N;
    n = 0:N-1;
    x = cos(0.4*pi*n) + A*cos(0.45*pi*n);
    X = zeros(size(omega));
    for k = 1:length(omega)
        X(k) = sum(x .* exp(-1j * omega(k) * n));
    end
    figure;
    plot(omega, abs(X),'LineWidth',1.2);
    xlabel('Frequency (rad/sample)'); ylabel('|X(e^{jω})|');
    title(sprintf('DTFT Magnitude, A=%.1f, N=%d', A, N));
    grid on;
end

%% Section 2.1: Synthetic Signals MRANC
% Parameters
N = 1000; fs = 500; t = (0:N-1)/fs;
% Desired signal
d = sin(2*pi*5*t + rand*2*pi);
% Interference
SIR_in = 0; A_int = 10^(-SIR_in/20);
i = A_int * sin(2*pi*50*t + rand*2*pi);
x = d + i;
% References
Z = [sin(2*pi*50*t); cos(2*pi*50*t)];
% Parameter combinations
mu_vals = [0.0005,0.005,0.01,0.01,0.01];
M_vals  = [0,      0,    5,   50,   550];
results = nan(length(mu_vals),5);

for k = 1:length(mu_vals)
    mu = mu_vals(k); M = M_vals(k);
    [s_out,~] = mranc(x,Z,M,mu);
    start_idx = round(0.7*N)+1;
    res_no = s_out(start_idx:end)-d(start_idx:end);
    SIR_out = 10*log10(var(d(start_idx:end))/var(res_no));
    orig_no = x(start_idx:end)-d(start_idx:end);
    SIR_in_act = 10*log10(var(d(start_idx:end))/var(orig_no));
    SIR_imp = SIR_out - SIR_in_act;
    err = abs(s_out-d); err_s = movmean(err,50);
    thresh = 1.5*mean(err_s(start_idx:end));
    idx_conv = find(err_s<thresh,1);
    conv_t = nan; if ~isempty(idx_conv), conv_t = idx_conv; end
    results(k,:) = [mu,M,SIR_out,SIR_imp,conv_t];
    % Plot for each
    figure;
    subplot(3,1,1); plot(t,d,'b',t,x,'r',t,s_out,'g','LineWidth',1.2);
    xlabel('Time (s)'); ylabel('Amplitude'); title(sprintf('MRANC, μ=%.4f, M=%d',mu,M)); legend('d','x','s'); grid on;
    subplot(3,1,2); plot(t,err,'r',t,err_s,'b','LineWidth',1.2); xlabel('Time (s)'); ylabel('Error'); grid on;
    if ~isnan(conv_t), hold on; plot([t(conv_t) t(conv_t)],ylim,'k--'); end
    subplot(3,1,3); [Px,f]=pwelch(x,hamming(256),128,1024,fs); [Ps,~]=pwelch(s_out,hamming(256),128,1024,fs);
    plot(f,10*log10(Px),'r',f,10*log10(Ps),'g','LineWidth',1.2); xlabel('Hz'); ylabel('dB/Hz'); legend('x','s'); grid on; xlim([0 100]);
end

% Display table
fprintf('--------------------------------------------------------------\n');
fprintf('| μ      | M    | SIR_out(dB) | SIR_imp(dB) | ConvTime |\n');
fprintf('--------------------------------------------------------------\n');
for k = 1:size(results,1)
    fprintf('| %.4f | %4d | %10.2f | %10.2f | %8.0f |\n', results(k,:));
end
fprintf('--------------------------------------------------------------\n');

%% MRANC Function
function [s, W] = mranc(x, Z, M, mu)
    [R,N] = size(Z); W = zeros(R*(M+1),1); s = zeros(1,N);
    ZZ = [];
    for r=1:R
        ZZ = [ZZ; toeplitz([Z(r,1) zeros(1,M)], Z(r,:))];
    end
    for n=M+1:N
        zvec = ZZ(:,n); y = W' * zvec;
        s(n) = x(n) - y;
        W = W + 2*mu*s(n)*zvec;
    end
end
