% %
% Zhiguo Wang and Bing Zhang
% 2013.06.3
% 
% Wang, Zhiguo, Bing Zhang, and Jinghuai Gao. 
% The residual phase estimation of a seismic wavelet using a R¨¦nyi divergence-based criterion. 
% Journal of Applied Geophysics 106 (2014): 96-105.
% %
function H=Renyientropyfun(x,L,a);
y=sort(x);
o1=ones(1,L)/L;
w=[-o1,o1];
o2=ones(1,2*L)/2/L;
u=conv(w,y);
v=(conv(o2,y));
v1 =normpdf(v,0,var(v));
J=(u.*v1).^(1-a);
H=real(log(sum(J)/length(J))/(a-1));


