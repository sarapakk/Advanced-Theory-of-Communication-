function [gsrvl ,gsrv2]=gngauss(m,sgma) 
if nargin == 0, m=O;
    sgma=1;
elseif nargin == 1 , sgma=m;
    m=0;
end;
u=rand;
z=sgma*(sqrt(2*log(1/(1-u))));
u=rand;
gsrvl=m+z*cos(2*pi*u);
gsrv2=m+z*sin(2*pi*u); 
