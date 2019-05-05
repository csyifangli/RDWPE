% %
% Zhiguo Wang and Bing Zhang
% 2013.06.3
% 
% Wang, Zhiguo, Bing Zhang, and Jinghuai Gao. 
% The residual phase estimation of a seismic wavelet using a R¨¦nyi divergence-based criterion. 
% Journal of Applied Geophysics 106 (2014): 96-105.
% %
function s=rickerfun(dt,L,f,fa);
t=-L*dt:dt:L*dt;
s0=(1-2*(pi*f)^2*t.^2).*exp(-(pi*f)^2*t.^2);
s00=hilbert(s0);
s000=real(s00)*cos(fa)+imag(s00)*sin(fa);
s=s000;



