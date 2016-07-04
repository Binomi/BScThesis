%% Simulation Ti:Sa femtosecond oscillator

N=2^13; %2^12; %number of time points
tmin=0; % start time
tmax=5e-12; % stop time (s)
t0=tmax/2; % center the pulse in the time window (s)
dt = (tmax - tmin)/(N-1); % time step
t=(tmin : dt : tmax); % Time scale

nz=5000; % round trips

% Parameters and constants
c=3e8; % (m/ps)
h=6.626e-34;

lambda = 800e-9; % Center wavelength [m]
f0=c/lambda; % center freq

% Fiber parameters
n=1.76;
n2=0.3*1e-4/1e15; %m^2/W
x=0.004; %m
beta=2*pi*n2*x/(lambda*n);% m^2/W (self-phase modulation)
gamma=0.24;   %modulation depth of sat.absorption
sigma=80;  % inverse saturaion intensity
alpha=0.074; % gain
rho=0.033; % loss output coupling

D0=-50e-30;

sigma_g=41e-24; %m^2
tf=2.5e-15;
ep=sigma_g*tf/(h*f0*beta);
Es=tf/(ep*beta);  %Gain Saturation

Tcav=12.8e-9;
Tr=3.2e-6;
alpha_m=0.33;
w0=10e-6; %m radius of pump-beam
Ip=4.6/(pi*w0^2*n); %Pump power normalized
v_a=c/532e-9;
sigma_a=20e-24;
P=Ip*sigma_a/(h*v_a)*Tcav;

% Start parameters
FW=25e-15;%3.74e-12; % full temporal width (s) at 1/e of maximum of Gaussian envelope
a=4/FW^2;
P0=1; %1e-3  %4; %130; % intensity (W/m^2), but normally peak power
A=sqrt(P0); % A is the amplitude of the wave to be propagated in the NLSE

u0=A*(exp(-a.*(t-t0).^2))+A*exp(-a*(t-t0-500e-15).^2+1i*pi);%.*exp(-1i*1e12*(t-t0));
U=u0;

df=1/tmax; %frequency increment
f=[-df*N/2:df:df*(N/2-1)];
f=fftshift(f);
wav = c./(f+c/lambda)*1e9; % wavelength scale of spectrum (nm)

% Dispersion
D=(-1i*D0 - tf^2)*(2*pi*f).^2;
disptot=exp(D);

ergebnis=zeros(N,nz);
alpha_hist=zeros(nz,1);
for k=1:nz
    E0=1/beta*sum(abs(U).^2)*dt; % pulse energy
    
    % Gain dynamics
    alpha = alpha * exp(-E0/Es - Tcav/Tr - P) + alpha_m*P/(P+E0/Es+Tcav/Tr)*(1-exp(-E0/Es-Tcav/Tr - P));
    alpha_hist(k)=alpha;
    alpha2=-1/beta*cumsum(abs(U).^2)/Es*dt*alpha;
    alpha2=alpha2-mean(alpha2);
    
    U=exp(alpha +alpha2 - rho ...
        -1i*abs(U).^2 ...
        - gamma./(1+sigma*abs(U).^2)).*U; % Nonlinearity
    U=fft(U);
    U=(disptot).*U; %dispersion step
    U=ifft(U);
    
    % save
    ergebnis(:,k)=U.';
end

figure
imagesc(1:nz,t,abs(ergebnis))

figure
imagesc(1:nz,fftshift(f),fftshift(abs(fft(ergebnis)),1))

