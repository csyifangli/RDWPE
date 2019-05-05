% %
% Zhiguo Wang and Bing Zhang
% 2013.12.5
% 
% Wang, Zhiguo, Bing Zhang, and Jinghuai Gao. 
% The residual phase estimation of a seismic wavelet using a R¨¦nyi divergence-based criterion. 
% Journal of Applied Geophysics 106 (2014): 96-105.
% %

clc; 
clear all; 
close all;
dt=0.001;
L=80;
for k=1:100;
a(k)=(-pi/2)+(k-1)*pi/100;
w(k,:)=rickerfun(dt,L,50,a(k)); 
P=ones(1,length(w(k,:)));
R1(k)=Renyientropyfun(w(k,:),10,0.1);
R2(k)=Renyientropyfun(w(k,:),10,0.3);
R3(k)=Renyientropyfun(w(k,:),10,0.5);
R4(k)=Renyientropyfun(w(k,:),10,0.7);
R5(k)=Renyientropyfun(w(k,:),10,0.9);

W(k,:)=abs(hilbert(w(k,:)));
Sim(k)=(sum(w(k,:).^1.*W(k,:))/(sum(w(k,:).^2)*sum(W(k,:).^2))^0.5);
K(k)=kurtosis(w(k,:));
G1(k)=sum(w(k,:).^2.*w(k,:))/(sum(w(k,:).^2)*sum(w(k,:).^2))^0.5;
G2(k)=sum(w(k,:).^2.*P)/(sum(w(k,:).^2)*sum(P))^0.5;
SK(k)=G1(k)^2/(G2(k)^2);

end

figure(1)
hold on
axis([-2,2,0,1]);
plot(a,R1/max(R1),'r')
plot(a,R2/max(R2),'g')
plot(a,R3/max(R3),'b')
plot(a,R4/max(R4),'y')
plot(a,R5/max(R5),'c')
plot(a,K/max(K),'k')
plot(a,Sim/max(Sim))
plot(a,SK/max(SK),'y')
xlabel('Phase (radian)')
ylabel('Normalized amplitude')
legend('R0.1','R0.3','R0.5','R0.7','R0.9','Kurtosis','Similarity','Skewness')
