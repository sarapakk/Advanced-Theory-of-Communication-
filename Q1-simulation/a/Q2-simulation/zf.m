clear all
clc
x = [1/12, 1/2, 5/6, 1/2, 1/12]; 
SNR = 0:12;
taps = 5; 
j = ftZF(x, taps, 0);

q = conv(j, x);
seqLength = 100000;
I = 2*randi([0 1],1,seqLength)-1;
n = ceil(length(q)/2);
while n+floor(length(q)/2) <= length(I)
    u(n-floor(length(q)/2)) = dot(q, I(n-floor(length(q)/2):n+floor(length(q)/2)));
    n = n + 1;
end
I = I(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);

temp = conv(x, j);
temp = conv(temp, j);  
Pe5 = zeros(size(SNR));
Pe_t5 = zeros(size(SNR));
for i = 1:length(SNR)
  
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    y = u + nl;
    I_ = 2*(real(y) > 0)-1;
    Pe5(i) = sum(I ~= I_)/length(I);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    for j = 0:2^(length(qq))-1
        Pe_t5(i) = Pe_t5(i) + 1/2^(length(qq))*qfunc((1 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    
    SNR(i);
end



plot(SNR, Pe5, 'r')
hold on
plot(SNR, Pe_t5, 'r.-')

fprintf('ZF with 9 taps\n');
clear u
taps = 9; 

j = ftZF(x, taps, 0);

q = conv(j, x);

seqLength = 100000;
I = 2*randi([0 1],1,seqLength)-1;
n = ceil(length(q)/2);
while n+floor(length(q)/2) <= length(I)
    u(n-floor(length(q)/2)) = dot(q, I(n-floor(length(q)/2):n+floor(length(q)/2)));
    n = n + 1;
end
I = I(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);

temp = conv(x, j);
temp = conv(temp, j); 

Pe9 = zeros(size(SNR));
Pe_t9 = zeros(size(SNR));
for i = 1:length(SNR)
    
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    y = u + nl;
    I_ = 2*(real(y) > 0)-1;
    Pe9(i) = sum(I ~= I_)/length(I);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    for j = 0:2^(length(qq))-1
        Pe_t9(i) = Pe_t9(i) + 1/2^(length(qq))*qfunc((1 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    
    SNR(i);
end


plot(SNR, Pe9, 'g')
plot(SNR, Pe_t9, 'g.-');

fprintf('ZF with 13 taps\n');
clear u
taps = 13; 
j = ftZF(x, taps, 0);
q = conv(j, x);


seqLength = 100000;
I = 2*randi([0 1],1,seqLength)-1;

n = ceil(length(q)/2);
while n+floor(length(q)/2) <= length(I)
    u(n-floor(length(q)/2)) = dot(q, I(n-floor(length(q)/2):n+floor(length(q)/2)));
    n = n + 1;
end
I = I(ceil(length(q)/2):ceil(length(q)/2)+length(u)-1);

temp = conv(x, j);
temp = conv(temp, j); 

Pe13 = zeros(size(SNR));
Pe_t13 = zeros(size(SNR));
for i = 1:length(SNR)
    
    snr = 10^(SNR(i)/10);
    N0 = 1/(2*snr);
    sigmaNoise = 2 * N0 * temp(ceil(length(temp)/2));
    nl = complex(sqrt(sigmaNoise/2)*randn(1, length(u)), sqrt(sigmaNoise/2)*randn(1, length(u)));
    y = u + nl;
    I_ = 2*(real(y) > 0)-1;
    Pe13(i) = sum(I ~= I_)/length(I);
    q0 = q(ceil(length(q)/2));
    qq = q(q ~= q0);
    for j = 0:2^(length(qq))-1
        Pe_t13(i) = Pe_t13(i) + 1/2^(length(qq))*qfunc((1 + dot((2*de2bi(j, length(qq))-1),qq) ) /sqrt(sigmaNoise/2));
    end
    
    SNR(i);
end


plot(SNR, Pe13, 'b')
plot(SNR, Pe_t13, 'b.-')
axis([0 12 0 0.5])

grid on
legend('5 taps sim', ' 5 taps theory',    ' 9 taps sim', ' 9 taps theory',   '13 taps sim', '13 taps theory')
xlabel('SNR-db')
ylabel('Pe')
