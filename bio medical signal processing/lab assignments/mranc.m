function [s, W] = mranc(x, Z, M, mu)
    % Multi-Reference Adaptive Noise Canceling (MRANC)
    % Inputs:
    %   x - Primary input signal (desired + noise)
    %   Z - Reference signals (matrix with multiple reference signals)
    %   M - Filter order
    %   mu - Step size parameter
    % Outputs:
    %   s - Output signal after noise cancellation
    %   W - Final adaptive filter weights
    
    % Initialize variables
    [R, N] = size(Z);
    W = zeros(R*(M+1), 1); % Column vector of filter weights
    s = zeros(1, N);       % Output signal
    
    % Generate Toeplitz matrix for reference signals
    ZZ = [];
    for r = 1:R
        ZZ = [ZZ; toeplitz([Z(r, 1), zeros(1, M)], Z(r, :))];
    end
    
    % Adaptive filtering loop
    for t = M+1:N
        % Extract current slice of reference signals
        z_vec = ZZ(:, t);
        
        % Compute filter output (estimated noise)
        y = W' * z_vec;
        
        % Compute error signal (desired signal estimate)
        s(t) = x(t) - y;
        
        % Update filter weights
        W = W + 2*mu * s(t) * z_vec;
    end
end