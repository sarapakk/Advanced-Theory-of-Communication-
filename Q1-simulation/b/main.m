clear all
clc
syms u x
SNRindB2 = 0:1:12;
%%%%%
a=sqrt((2-sqrt(2))*SNRindB2);
b=sqrt((2+sqrt(2))*SNRindB2);
l=a.*b;
io=(exp(u))./sqrt(2*pi.*l);
u=a.*x;
ss=x.*exp(-(x.^2+a.^2)./2).*io;
for j=1:length(ss)
q=int(ss(j),b(j),inf);
end
u=a.*b;
pb=q-((1/2)*io.*exp(-(x.^2+a.^2)./2));

SNRindBl = 0:1:12; 
for i=1 :length(SNRindBl),
    srnld_err_prb(i)=cm_sm34(SNRindBl(i)); 

end; 
    
    for i=1 :length(SNRindB2), 
        SNR=exp(SNRindB2(i)*log(10)/10);
        theo_err_prb(i)=2* qfunc( sqrt(SNR));
        
    end;
    figure
        semilogy(SNRindBl, srnld_err_prb ,'mx-','LineWidth',2) ;
        hold on
        semilogy(SNRindB2,theo_err_prb,'bs-','Linewidth',2);
        
        axis([0 12 10^-3 1])
grid on
xlabel('Eb/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for dqpsk modulation')
