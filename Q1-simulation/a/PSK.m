% Bit Error Rate for 16-PSK modulation

clear all
clc
N = 10^5; % number of symbols
M = 16; k = log2(M); 
thetaMpsk = [0:M-1]*2*pi/M; % reference phase values

Eb_N0_dB = [0:12]; % multiple Es/N0 values
Es_N0_dB = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = [0:M-1];
map = bitxor(ref,floor(ref/2));
[tt ind] = sort(map);
ipPhaseHat = zeros(1,N);
for ii = 1:length(Eb_N0_dB)

 % symbol generation
 ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
 bin2DecMatrix = ones(N,1)*(2.^[(k-1):-1:0]) ; % conversion from binary to decimal
 ipBitReshape = reshape(ipBit,k,N).'; % grouping to N symbols having k bits each
 ipGray = [sum(ipBitReshape.*bin2DecMatrix,2)].'; % decimal to binary
 
 % Gray coded constellation mapping
 ipDec = ind(ipGray+1)-1; % bit group to constellation point 
 ipPhase = ipDec*2*pi/M; % conversion to phase 
 ip = exp(j*ipPhase);% modulation 
 s = ip; 
 
 % noise

 n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance 
 
 y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise
 % finding the phase from [-pi to +pi]
 opPhase = angle(y); 
 
 opPhase(find(opPhase<0)) = opPhase(find(opPhase<0)) + 2*pi;
 % rounding the received phase to the closest constellation
 ipPhaseHat = 2*pi/M*round(opPhase/(2*pi/M)) ;
 
 % this is to enable comparison with the transmitted phase
 ipPhaseHat(find(ipPhaseHat==2*pi)) = 0;
 ipDecHat = round(ipPhaseHat*M/(2*pi));

 % Decimal to Gray code conversion
 ipGrayHat = map(ipDecHat+1); 
 ipBinHat = dec2bin(ipGrayHat,k) ;

 % converting binary string to number
 ipBinHat = ipBinHat.';
 ipBinHat = ipBinHat(1:end).';
 ipBinHat = str2num(ipBinHat).' ;
 
 % counting errors
 nBitErr(ii) = size(find([ipBit- ipBinHat]),2); % couting the number of errors

end 
simBer = nBitErr/(N*k);
theoryBer = (1/k)*erfc(sqrt(k*10.^(Eb_N0_dB/10))*sin(pi/M));

close all;
figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2)
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2)
axis([0 12 10^-3 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-PSK modulation')