clear all
clc
SNR = 0:12;
x = [1/12, 1/2, 5/6, 1/2, 1/12];k1 = 2;k2 = 2;
sl = 100000;
I = 2*randi([0 1], 1, sl)-1;
R = zeros(length(-k1-k2:k1+k2), 1);
R(1:length(x)) = x;
R = circshift(R, [(length(R)-length(x))/2, 0]);
R_IYmat = ftdfe(R, k1, k2);
Pe5 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    [R_Ymat, R_Y] = ft2dfe (x, k1, snr);
    R_Y = R_Y.';
    A = [R_Ymat, R_IYmat; R_IYmat', eye(k2)];
    Ans = zeros(length(-k1:k2), 1);
    Ans(1:length(1+k2:k1+k2+1)) = R(1+k2:k1+k2+1);
   
    d = A^(-1)*Ans;
    
    FF = d(1:k1+1);
    FB = d(k1+2:end);

    q = conv(FF, x);
    q_index = -k1-floor(length(x)/2): -k1-floor(length(x)/2) + length(q)-1;
    q0 = q(q_index == 0);
    qq = q(q ~= q0);
    
    I_ = zeros(size(I));
    I_(1:k2) = I(1:k2);
    n = k2+1;
    sigmaNoise = 2*N0*x(ceil(length(x)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(I)), sqrt(sigmaNoise/2)*randn(1, length(I)));
    nl = filter(FF, 1, nl);
    while n - q_index(1) <= length(I)
        temp = dot(q, I(n-q_index)) + nl(n);
        y = temp + dot(I_(n-1:-1:n-k2), FB);
        I_(n) = 2*(real(y) > 0) - 1;
        n = n + 1;
    end
    I__ = I(1:length(I_));

    Pe5(i) = sum(I__ ~= I_)/length(I);
    SNR(i);
end



semilogy(SNR, Pe5, 'r')
hold on

k1 = 4;
k2 = 4;

sl = 100000;
I = 2*randi([0 1], 1, sl)-1;

R = zeros(length(-k1-k2:k1+k2), 1);
R(1:length(x)) = x;
R = circshift(R, [(length(R)-length(x))/2, 0]);

R_IYmat = ftdfe(R, k1, k2);


Pe9 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    
   
    [R_Ymat, R_Y] = ft2dfe (x, k1, snr);
    R_Y = R_Y.';
    A = [R_Ymat, R_IYmat; R_IYmat', eye(k2)];
    Ans = zeros(length(-k1:k2), 1);
    Ans(1:length(1+k2:k1+k2+1)) = R(1+k2:k1+k2+1);
    
    
    d = A^(-1)*Ans;
    
   
    FF = d(1:k1+1);
    FB = d(k1+2:end);

    q = conv(FF, x);
    q_index = -k1-floor(length(x)/2): -k1-floor(length(x)/2) + length(q)-1;
    q0 = q(q_index == 0);
    qq = q(q ~= q0);
    
 
    I_ = zeros(size(I));
    I_(1:k2) = I(1:k2);
    n = k2+1;
    sigmaNoise = 2*N0*x(ceil(length(x)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(I)), sqrt(sigmaNoise/2)*randn(1, length(I)));
    nl = filter(FF, 1, nl);
    while n - q_index(1) <= length(I)
        temp = dot(q, I(n-q_index)) + nl(n);
        y = temp + dot(I_(n-1:-1:n-k2), FB);
        I_(n) = 2*(real(y) > 0) - 1;
        n = n + 1;
    end
    I__ = I(1:length(I_));

    Pe9(i) = sum(I__ ~= I_)/length(I);
    SNR(i);
end

semilogy(SNR, Pe9, 'g')
hold on


k1 = 6;
k2 = 6;

% generating the BPSK sequence
sl = 100000;
I = 2*randi([0 1], 1, sl)-1;

% cross-correlation matrix
R = zeros(length(-k1-k2:k1+k2), 1);
R(1:length(x)) = x;
R = circshift(R, [(length(R)-length(x))/2, 0]);

R_IYmat = ftdfe(R, k1, k2);

Pe13 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    
    % auto-correlation matrix
    [R_Ymat, R_Y] = ft2dfe (x, k1, snr);
    R_Y = R_Y.';
    A = [R_Ymat, R_IYmat; R_IYmat', eye(k2)];
    Ans = zeros(length(-k1:k2), 1);
    Ans(1:length(1+k2:k1+k2+1)) = R(1+k2:k1+k2+1);
    
    
    d = A^(-1)*Ans;
    
    % feedback and feedforward coefficients
    FF = d(1:k1+1);
    FB = d(k1+2:end);

    q = conv(FF, x);
    q_index = -k1-floor(length(x)/2): -k1-floor(length(x)/2) + length(q)-1;
    q0 = q(q_index == 0);
    qq = q(q ~= q0);
    
    
    I_ = zeros(size(I));
    I_(1:k2) = I(1:k2);
    n = k2+1;
    sigmaNoise = 2*N0*x(ceil(length(x)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(I)), sqrt(sigmaNoise/2)*randn(1, length(I)));
    nl = filter(FF, 1, nl);
    while n - q_index(1) <= length(I)
        temp = dot(q, I(n-q_index)) + nl(n);
        y = temp + dot(I_(n-1:-1:n-k2), FB);
        I_(n) = 2*(real(y) > 0) - 1;
        n = n + 1;
    end
    I__ = I(1:length(I_));

    Pe13(i) = sum(I__ ~= I_)/length(I);
    SNR(i);
end


semilogy(SNR, Pe13, 'b')

axis('square')
grid on
legend(' 5 tap', '9', ' 13')
xlabel('SNR-db')
