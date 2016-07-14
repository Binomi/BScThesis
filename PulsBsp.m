a1=1;
a2=sqrt(1);
tau=10;
phi=0.5*pi;
dt=0.1;
t=(-100:dt:100)';
N=length(t);
df=1/(max(t)-min(t));
f=df*(-(N-1)/2:(N-1)/2)';
y=a1*exp(-(1.1779*t).^2)+a2*exp(-(1.1779*(t-tau)).^2+1i*phi);
spectrum=fftshift(abs(fft(y)).^2);
spectrum=spectrum./max(spectrum);
spectrum=circshift(spectrum,[0,0]);
autocorr=fftshift(fft(spectrum));
autocorr=autocorr./max(abs(autocorr));
% autocorr2=zeros(N,1);
% for ii=1:N
%     autocorr2(ii)=conj(y)*circshift(y,[0,ii-1])';
% end
% autocorr2=fftshift(autocorr2);
% autocorr2=autocorr2./max(autocorr2);
figure
plot(t,abs(y).^2)
xlim([-15,15])

figure
plot(f,spectrum)
xlim(1.2*[-1,1])

t2=t(1:2:end);
phase=angle(autocorr(1:2:end));
phiTrend=polyfit(t2((t2<2) & (t2>-2)),phase((t2<2) & (t2>-2)),1);
phase2=phase-phiTrend(1)*t2;

figure
hold
plot(t,abs(autocorr),t2,phase2)
%plot(x,abs(autocorr2),x,angle(autocorr2))
xlim([-15,15])

[~,loc]=findpeaks(abs(autocorr(1:2:end)),'minPeakProminence',0.1);
phi1=phase2(loc);
figure
plot(phi1)