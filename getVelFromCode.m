%% This script should be used after readFileWithSetTime to get the velocity
% of a meteor when you know the coordinates of the meteor in the raw data
% set needed constants
freq = 450e6;       %radar freq    (Hz)
c = 2.99792458e8;   %speed of light(m/s)
wavelen = c/freq;   %wavelength in (m)
dtPixel = 5e-7;     %pixel size in time (s)
velRange = 7;      %velocity range to check (+/-) (km/s)
velGuess = -1.5;       %velocity guess (based on RTI fit)
nVels = 1000;        %number of velocities to try

wideRawRx = rawdat(:,2:2:end,irec);
narrRawRx = rawdat(:,1:2:end,irec);
wideTx = rawtx(:,2:2:end,irec);
narrTx = rawtx(:,1:2:end,irec);

%plot the wide beam
imagesc(db(abs(wideRawRx)))
xlabel('x pixel (wideRaw(y,x))')
ylabel('y pixel (wideRaw(y,x))')
title('Wide Beam Raw Data (dB)')


%now the user should zoom in on the meteor
%get the approximate region where the meteor is
xStart = 306; %600;  %CHANGE THIS BASED ON THE COORDINATES
yStart = 504; %450;  %CHANGE THIS BASED ON THE COORDINATES
xEnd = 409; %628;    %CHANGE THIS BASED ON THE COORDINATES
yEnd = 640;  %;    %CHANGE THIS BASED ON THE COORDINATES
if (yEnd-yStart) < size(wideTx,1)
    yEnd = yStart + size(wideTx,1)+1;
end
nPulses = xEnd-xStart+1;
%now do a convolution with the tx signal for each pulse
convMagsAll = zeros(nPulses,nVels);
velsForEachPulse = zeros(nPulses,1);
for pulseNumMet = 1:nPulses
    %get the actual pulse index
    pulseIndex = pulseNumMet+xStart-1;
    %isolate the pulse we care about
    metPulse = wideRawRx(yStart:yEnd,pulseIndex);
    metPulse = flipud(metPulse);
    txPulse = wideTx(:,pulseIndex)';
    %generate possible velocities (m/sec)
    velPossible = linspace(velGuess-velRange,velGuess+velRange,nVels)*1000;
    convMags = zeros(1,numel(velPossible));
    convInds = zeros(1,numel(velPossible));
    %do a convolution with the code
    for velIndex = 1:numel(velPossible)
        velTest = velPossible(velIndex);
        %get the estimated delta phase for each pixel
        dPhases = linspace(0,numel(txPulse)-1,numel(txPulse))*...
            velTest*dtPixel/wavelen*2*pi;
        txPulseMod = txPulse.*exp(1i*dPhases);
        [convMags(velIndex),convInds(velIndex)]= max(abs(conv(metPulse,txPulseMod,'same')));
    end
    %put the magnitudes into the convolution array
    convMagsAll(pulseNumMet,:) = convMags;
    %get the velocities
    [~,velMaxIndex] = max(convMags);
    velsForEachPulse(pulseNumMet) = velPossible(velMaxIndex);
end

%plot the data
figure()
plot(velsForEachPulse)
xlabel('Pulse Number')
ylabel('Velocity (m/s)')

%make a nicer plot of the raw data
backIndeces = 9;
figure()
imagesc(tms(1,xStart-backIndeces:end)-tms(1,xStart-backIndeces),rng/1000, ...
    db(abs(wideRawRx(:,xStart-backIndeces:end))/median(abs(wideRawRx(:)))))
% imagesc(tms(1,xStart-backIndeces:end)-tms(1,xStart-backIndeces),rng/1000, ...
%     db(abs(narrRawRx(:,xStart-backIndeces:end))/median(abs(narrRawRx(:)))))
set(gca,'ydir','normal')
set(gca,'ylim',[98,130])
set(gca,'xlim',[0,0.36])
fontsize = 20;
set(gca,'fontsize',fontsize)
set(gca,'clim',[-3,20])
colorbar('fontsize',fontsize)
xlabel('Time (sec)','fontsize',fontsize)
ylabel('Range (km)','fontsize',fontsize)
title('Raw Data')
figure()
plot(velPossible,convMagsAll(40,:),'linewidth',2)
set(gca,'fontsize',fontsize)
