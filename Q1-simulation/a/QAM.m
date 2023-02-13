% symbol error rate for 16-QAM modulation
clear all
clc
N = 2*10^5; 
alpha16qam = [-3 -1 1 3]; k=4;
Eb_N0_dB = [0:12]; % multiple Es/N0 values
Es_N0_dB = Eb_N0_dB + 10*log10(k); 
ipHat = zeros(1,N);
for ii = 1:length(Es_N0_dB)
 ip = randsrc(1,N,alpha16qam) + j*randsrc(1,N,alpha16qam);
 s = (1/sqrt(10))*ip; % normalization of energy to 1
 n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

 y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

 % demodulation
 y_re = real(y); % real part
 y_im = imag(y); % imaginary part

 ipHat_re(find(y_re< -2/sqrt(10))) = -3;
 ipHat_re(find(y_re > 2/sqrt(10))) = 3;
 ipHat_re(find(y_re>-2/sqrt(10) & y_re<=0)) = -1;
 ipHat_re(find(y_re>0 & y_re<=2/sqrt(10))) = 1;
 ipHat_im(find(y_im< -2/sqrt(10))) = -3;
 ipHat_im(find(y_im > 2/sqrt(10))) = 3;
 ipHat_im(find(y_im>-2/sqrt(10) & y_im<=0)) = -1;
 ipHat_im(find(y_im>0 & y_im<=2/sqrt(10))) = 1;
 ipHat = ipHat_re + j*ipHat_im;
 nErr(ii) = size(find([ip- ipHat]),2);
 
end

simBer = nErr/N;
theoryBer = 3/2*erfc(sqrt(0.1*(10.^(Es_N0_dB/10))));
close all
figure
semilogy(Eb_N0_dB,theoryBer,'b.-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','Linewidth',2);
axis([0 12 10^-3 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16-QAM modulation') 