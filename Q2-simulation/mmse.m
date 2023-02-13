clc
clear all
x = [1/12, 1/2, 5/6, 1/2, 1/12]; SNR = 0:12;
sl = 100000;
t= 2*randi([0 1],1,sl)-1;
num = 5;
k = zeros(num, 1);
k(1:length(x)) = x;
k = circshift(k, [(num-length(x))/2, 0]);
Pe5 = zeros(size(SNR));
Pe_th5 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    R_Ymat = ft (x, num, snr);
    d = R_Ymat^(-1)*k;
    q = conv(d, x);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    n = ceil(length(q)/2);
    while n+floor(length(q)/2) <= length(t)
        u(n-floor(length(q)/2)) = dot(q, t(n-floor(length(q)/2):n+floor(length(q)/2)));
        n = n + 1;
    end
    t = t(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);
    
    temp = conv(x, d);
    temp = conv(temp, d);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    
    y = u + nl;
    for j = 0:2^(length(qq))-1
        Pe_th5(i) = Pe_th5(i) + 1/2^(length(qq))*qfunc((q0 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    clear u
    
    I_ = 2*(real(y) > 0)-1;
    Pe5(i) = sum(t ~= I_)/length(t);
    SNR(i);
end

semilogy(SNR, Pe5, 'r')
hold on
semilogy(SNR, Pe_th5, 'r.-');


sl = 100000;
t = 2*randi([0 1],1,sl)-1;

num = 9; 

k = zeros(num, 1);
k(1:length(x)) = x;
k = circshift(k, [(num-length(x))/2, 0]);
Pe9 = zeros(size(SNR));
Pe_t9 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    
    R_Ymat = ft(x, num, snr);
    
    d = R_Ymat^(-1)*k;
    q = conv(d, x);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    
    n = ceil(length(q)/2);
    while n+floor(length(q)/2) <= length(t)
        u(n-floor(length(q)/2)) = dot(q, t(n-floor(length(q)/2):n+floor(length(q)/2)));
        n = n + 1;
    end
    t = t(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);

    temp = conv(x, d);
    temp = conv(temp, d);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    
    y = u + nl;

    for j = 0:2^(length(qq))-1
        Pe_t9(i) = Pe_t9(i) + 1/2^(length(qq))*qfunc((q0 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    clear u
    
    I_ = 2*(real(y) > 0)-1;
    Pe9(i) = sum(t ~= I_)/length(t);
    SNR(i);
end


semilogy(SNR, Pe9, 'g')
semilogy(SNR, Pe_t9, 'g.-');
sl = 100000;
t = 2*randi([0 1],1,sl)-1;

num = 13; 
k = zeros(num, 1);
k(1:length(x)) = x;
k = circshift(k, [(num-length(x))/2, 0]);
Pe13 = zeros(size(SNR));
Pe_t13 = zeros(size(SNR));
for i = 1:length(SNR)
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
   
    R_Ymat = ft (x, num, snr);
    
    d = R_Ymat^(-1)*k;
    q = conv(d, x);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    n = ceil(length(q)/2);
    while n+floor(length(q)/2) <= length(t)
        u(n-floor(length(q)/2)) = dot(q, t(n-floor(length(q)/2):n+floor(length(q)/2)));
        n = n + 1;
    end
    t = t(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);

    temp = conv(x, d);
    temp = conv(temp, d);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    
    y = u + nl;
    for j = 0:2^(length(qq))-1
        Pe_t13(i) = Pe_t13(i) + 1/2^(length(qq))*qfunc((q0 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    clear u
    
    I_ = 2*(real(y) > 0)-1;
    Pe13(i) = sum(t ~= I_)/length(t);
    SNR(i);
end

semilogy(SNR, Pe13, 'b')
semilogy(SNR, Pe_t13, 'b.-');
axis([0 12  0 0.5])
grid on
legend(' 5 taps sim', ' 5 taps '    ,'9 taps sim', '9 taps theory'    ,'13 taps sim', ' 13 taps theory')
xlabel('SNR-db')
ylabel('Pe')
