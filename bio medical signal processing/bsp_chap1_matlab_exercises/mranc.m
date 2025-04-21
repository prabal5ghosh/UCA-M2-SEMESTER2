function[s, W] = mracnc(x, Z, M, mu)
    % Initialize variables
    [R, N] = size(Z);
    W = zeros( R*(M + 1), N);
    s = zeros(1, N);



    % Generate Toeplitz matrix ZZ
    ZZ = [];
    for r = 1:R
        ZZ = [ZZ; toeplitz([Z(r, 1), zeros(1, M)], Z(r, :))];
    end

    for t = 1:N-1
        y(t) = W(:)' *  ZZ(:, t);
        s(t) = x(t) - y(t);
        W = W + 2*mu * s(t) * ZZ(:, t);
    end


























































