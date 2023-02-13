function d = ZF (x, taps, symParam)
% you can use symParam = 0 for symmetric ZF
    k = (taps-1)/2;
    L = (length(x)-1)/2;
    q_length = 2*(k+L)+1;
    
    x_ = x(end:-1:1);
    X = toeplitz([x_, zeros(1, 2*k)], zeros(1, taps));
    
    q = zeros(q_length, 1);
    q((q_length+1)/2) = 1;
    
    start = (q_length+1)/2-k-symParam;
    if start <= 0
        warning('symParam should be <= %d', (q_length+1)/2-k-1);
        start = 1;
    elseif start > (q_length+1)/2
        warning('symParam should be > %d', -k);
        start = (q_length+1)/2;
    end
    
    X = X(start:start+taps-1, :);
    q = q(start:start+taps-1);
    d = X^(-1) * q;
end