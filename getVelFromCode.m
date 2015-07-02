%% This script should be used after readFileWithSetTime to get the velocity
% of a meteor when you know the coordinates of the meteor in the raw data
% set needed constants
freq = 450e6;       %radar freq    (Hz)
c = 2.99792458e8;   %speed of light(m/s)
wavelen = c/freq;   %wavelength in (m)
dtPixel = 5e-7;     %pixel size in time (s)

wideRawRx = rawdat(:,2:2:end,irec);
narrRawRx = rawdat(:,1:2:end,irec);
wideTx = rawtx(:,2:2:end,irec);
narrTx = rawtx(:,1:2:end,irec);

%plot the wide beam
imagesc(db(abs(wideRawRx)))
xlabel('x pixel (wideRaw(y,x))')
ylabel('y pixel (narrRaw(y,x))')
title('Wide Beam Raw Data (dB)')

%now the user should zoom in on the meteor
%get the approximate region where the meteor is
xStart = 600;  %CHANGE THIS BASED ON THE COORDINATES
xEnd = 628;    %CHANGE THIS BASED ON THE COORDINATES
yStart = 450;  %CHANGE THIS BASED ON THE COORDINATES
yEnd = 610;    %CHANGE THIS BASED ON THE COORDINATES
if (yEnd-yStart) < size(wideTx,1)
    yEnd = yStart + size(wideTx,1)+1;
end
nPulses = xEnd-xStart+1;
%now do a convolution with the tx signal for each pulse
for pulseNumMet = 1:nPulses
    %get the actual pulse index
    pulseIndex = pulseNumMet+xStart-1;
    %isolate the pulse we care about
    metPulse = wideRawRx(yStart:yEnd,pulseIndex);
    txPulse = wideTx(:,pulseIndex)';
    %generate possible velocities
    velPossible = linspace(0,140,300)*1000;
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
    
    
    
end