clc;
clear all;
close all;
nsub=16;
cp=8;
nfft=64;
M=16;
%=========================================================
N=nfft+cp/2;
F=dftmtx(N)/sqrt(N);
%===========================================================
%============= Iterative =====================================
Ebno=0:2:22;
cl={'rx-','kx-','mx-'};bb{1}=1;
for it=1:2
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
h=rand(1,size(xf,2))+i*rand(1,size(xf,2));
Yp=awgn(diag(F*h(1,:)')*xf(1,:)',Ebno(nn));
Yn=awgn(diag(F*h(1,:)')*xf(2,:)',Ebno(nn));
%====== end of transmitter ==============
if it==1
yp=inv(diag(F*h(1,:)'))*Yp;
yn=inv(diag(F*h(1,:)'))*Yn;
yd=yp(cp/2+1:end)-(-yn(cp/2+1:end));
yd2=yp-(-yn);
yf2=(real(yd2)-i*imag(yd2));
yf=(real(yd)-i*imag(yd));
xdm=qamdemod(fft(yf),M);
Xf{it}=yf2;
else
    sz=sign(diag(Xf{it-1}));
    G=0.5*[1+F*sz*F';F*sz*F'-1];
    Tz=(G'*G).^-1*G';
  %==================================  
    yp=(.5*(1+F*sz*F'))*Yp;
    yn=(0.5*(F*sz*F'-1))*Yn;
    yd2=yp-(-yn);yf2=(real(yd2)-i*imag(yd2));
    Xf{it}=yf2;
    %==========================
    yd=yp(cp/2+1:end)-(-yn(cp/2+1:end));
    yf=(real(yd)-i*imag(yd));
    xdm=qamdemod(fft(yf),M);                
end     
[nr br(ii)]=symerr(xdm,X);
end
brf(nn)=mean(br)/nn;  
bri(nn)=berfading(Ebno(nn),'qam',M,1)/nn;
end
brf=brf.*(bb{1});
bb{1}=brf;
semilogy(Ebno,smooth(brf),cl{it},'Linewidth',2);hold on
end
xlim([0 20]);
grid on;
xlabel('---Ebno');
%ylim([10^-4 10^0]);
ylabel('--BER')
semilogy(Ebno,bri,'mx-','linewidth',2)
legend('Conven','Iterative','Lower bound',3);