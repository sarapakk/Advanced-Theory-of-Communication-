function [p]=cm_sm34(snr_in_dB)
 N=10000; E=1;  snr=10^(snr_in_dB/10); 
 sgma=sqrt(E/(4*snr));
 for i=1 :2*N, 
     temp=rand; 
     if (temp<0.5), dsource(i)=0; 
     else dsource(i)=1 ; 
     end;
 end;
 mapping=[0 1 3 2]; M=4;
 [diff_enc_output] = cm_dpske(E,M,mapping,dsource);
 for i=1 :N,
     [n(1) n(2)]=gngauss(sgma); 
     r(i,:)=diff_enc_output(i,:)+n; end;
 numoferr=0; 
 prev_theta=0; 
 for i=1 :N, theta=angle(r(i, 1 )+j*r(i,2));
     delta_theta=mod(theta-prev_theta,2*pi);
  if((delta_theta<pi/4)|(delta_theta> 7*pi/4)), 
      decis=[0 0];
      elseif (delta_theta<3*pi/4),
          decis=[0 1];
  elseif (delta_theta<5*pi/4) decis=[1 1 ];
  else decis=[1 0]; 
  end;
  prev_theta=theta;
  if ((decis(1)~=dsource(2*i-1 ))|(decis(2)~=dsource(2*i))),
      numoferr=numoferr+ 1 ; 
  end;
  end; p=numoferr/N;
  
  
      