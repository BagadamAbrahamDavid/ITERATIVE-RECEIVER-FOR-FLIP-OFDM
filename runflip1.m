clc;
clear all;
close all;
nsub=16;
cp=32;
nfft=64;
M=16;
%=================================================
N=nfft+cp/2;
F=dftmtx(N)/sqrt(N);
%=================================================
Ebno=0:2:20;
for nn=1:length(Ebno)
for ii=1:nsub
X=randi([0 M-1],nfft,1);
Xm=qammod(X,M);
x=ifft(Xm);
%================================================
xout=[x];
for l=1:length(xout)
    if xout(l)<0
        xp(l)=0;
        xn(l)=xout(l);
    else
        xp(l)=xout(l);
        xn(l)=0;
    end
end
xf=[[xp(nfft-cp/2+1:nfft) xp]; [xn(nfft-cp/2+1:nfft) xn]]; % ref:9 Fig:3 
rand('state',2)
nt=rand(size(xf))+i*rand(size(xf));
h=rand(size(xf))+i*rand(size(xf));
No=10^-(Ebno(nn)/10);
Yp=diag(F*h(1,:)')*xf(1,:)'+(No/2)*nt(1,:)';
Yn=diag(F*h(2,:)')*xf(2,:)'+(No/2)*nt(2,:)';
%====== end of transmitter ==============
yp=inv(diag(F*h(1,:)'))*Yp;
yn=inv(diag(F*h(2,:)'))*Yn;
yd=yp(cp/2+1:end)-(-yn(cp/2+1:end));
yf=(real(yd)-i*imag(yd));
xdm=qamdemod(fft(yf),M);
[nr br(ii)]=symerr(xdm,X);
end
brf(nn)=mean(br)/nn;  
end
semilogy(Ebno,smooth(brf),'x-');hold on
xlim([0 Ebno(end)]);
grid on;
xlabel('---Ebno');
ylabel('--BER')
