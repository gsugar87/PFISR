% calculate response in camera to star magnitudes

clear all; close all;

% camera specs: sCMOS, EMCCD, ICCD, perfect
framerate = [30 30 30 30 30];
lensd = 25e-3;
lensf = 1.4;
pixd = [6.5 6.5 16 10 10]*1e-6;     % microns
lambda = 600e-9;
QE = [0.72 0.72 0.95 0.5 1];
gain = [1 1 100 100 1];
readnoise = [1.2 1.2 1.0 0 0];
darknoise = [0.14 0.14 0.0001 0 0];
enf = [1.0 1.0 1.4 1.6 1.0]; % excess noise factor for EMCCD

background = 5000; % 5 kR at 600 nm (roughly)

% convert star magnitudes to brightness (flux)
M_sun = -26.77;   % Absolute bolo mag
F_sun = 1360; % W/m2
M = 5:0.1:18;

Flux = F_sun * 2.512.^(M_sun - M);

% assume we have an f/1, 50 mm lens. Area is then pi*(0.05/2)^2.

power = Flux * pi * (lensd/2)^2;

% convert watts to photons / sec, assuming 600 nm wavelength

h = 6.626e-34;
c = 3e8;
energyperphoton = h*c/lambda;

photonspersec = power / energyperphoton;

h1 = figure(1);
ax = axes;
set(ax,'fontsize',14);

colors = [{'k'} {'k--'} {'b'} {'r'} {'k:'}];

for k = 1:length(framerate),
    
    
    inttime = 1/framerate(k);
    photons = photonspersec * inttime;
    sig_electrons = photons * QE(k) * gain(k);
    
    % photon background
    aomega = (pi/4) * (pixd(k)*100/lensf)^2;
    bg_photons = background * 1e6/4/pi * aomega * inttime;
    
    if k == 2,
        sig_electrons = sig_electrons * 4;
        bg_photons = bg_photons * 4;
        readnoise(k) = readnoise(k) * 4;
        darknoise(k) = darknoise(k) * 4;
    end
    
    bg_electrons = QE(k) * bg_photons * gain(k);
    
    noise = sqrt( (enf(k)*gain(k))^2 * (sig_electrons/gain(k) + darknoise(k)^2) + readnoise(k)^2 );
    
    SNR = sig_electrons ./ noise;
    
    
    hand(k) = semilogy(ax,M-3,SNR,colors{k},'linewidth',1.2);
    hold(ax,'on');
    plot(ax,[M(1) M(end)]-3,[2 2],'r:');
    xlabel('Star Magnitude','Fontsize',16);
    ylabel('Signal to Noise Ratio','Fontsize',16);
    
    ind = find(SNR < 2,1,'first');
    plot(ax,[M(ind) M(ind)]-3,[min(SNR) max(SNR)],'b:');
    text(4,3,'SNR = 2','parent',ax,'fontsize',16);
    text(M(ind)+0.2-3,10+10*k,sprintf('%.1f Magnitude',M(ind)-3),'parent',ax,'fontsize',16);
    
    i = find(M>12,1,'first');
    fprintf('SNR at magnitude 10 = %.2f\n',SNR(i));
    
end

set(ax,'ylim',[1e-1 1e2]);
set(ax,'xlim',[2 15]);

%% for comparison, plot Hawkes equation 5.4.3

hSNR = logspace(-2,2,200);
R = 1e5;    % distance in meters
tau = 1/framerate(1);    % integration time
D = lensd;
n = hSNR.^2;     % number of signal photons
gamma = 0.9;    % lens transmission coeff.
eps = QE(1);

me = 49.9 - 2.5 * log10(n*R^2/(gamma*eps*D^2*tau));

hand(k+1) = plot(ax,me,hSNR,'m');

legend(hand,'sCMOS','sCMOS binned','EMCCD','ICCD','Perfect','Hawkes');

