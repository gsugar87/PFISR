%run the anayzeSacedMeteors script
analyzeSavedMeteors

%get constants
powerTx = 1.4e6; %Power in MegaWatts
ipp = 0.003; %inter pulse period in sec
returnDir = cd;
%get constants
cd('C:\Users\Glenn\Documents\MATLAB\Meteors\')
constants
%get back to original dir
cd(returnDir)

%% GET THE NORMALIZED RADAR POWER %%
%first we need to expand the "good indeces" of the radar for the meteor
% meteor 9 (startTime = 3.93 endTime = 4.23)
%1: 4.6, 4.8
%2: 3.32, 3.6
%3: 1.7, 2
%5: 2.35, 2.95
%6: 1.3, 1.9
%7: 1.9, 2.05
%9: 3.93, 4.23
%13: 0.84, 1.33
%18: 1.3, 2.65
startTimes = [4.63, 3.35, 1.728, 2.415, 1.375, 1.93, 3.82, 0.885,  1.365, 3.383,  0,   2.85];
endTimes =   [4.761, 3.57, 1.98, 2.916, 1.862, 2.04, 4.35, 1.262,  2.58, 3.686,  0.23, 3.02];
meteorIdentities = [1, 2,     3,     5,     6,    7,    9,    13,  18,  20, 24,    27]; %24
beamTypes =        [1, 1,     1,     1,     1,    1,    1,     1,   1,    2,  1,   1];  
extraTimeShifts = [0.2,-0.08,-0.08,-0.08, -0.02,-0.04,0.06,-0.1, 0.06,0.22,0.17,   0.18];
extraTimeShift = extraTimeShifts(meteorIdentities == meteorNumber); %+0.2; %-.08;
beamType = beamTypes(meteorIdentities == meteorNumber); %assume wide beam
startTime = startTimes(meteorIdentities == meteorNumber);
endTime = endTimes(meteorIdentities == meteorNumber);
%startTime = 0;
%endTime = 0.23;
if beamType == 1
    allPower = wideBeamPower;
    allTime = wideBeamTime;
    allRange = wideBeamRange;
else
    allPower = narrowBeamPower;
    allTime = narrowBeamTime;
    allRange = narrowBeamRange;
end
[~,radarStartIndexTime] = min(abs((allTime-allTime(1))*matlabTmsInSec - startTime));
[~,radarEndIndexTime] = min(abs((allTime-allTime(1))*matlabTmsInSec - endTime));
%get the time indeces that the meteor exists
radarTimeIndeces = linspace(radarStartIndexTime,...
    radarEndIndexTime, radarEndIndexTime-radarStartIndexTime + 1);

%now that we have the time indeces we want to plot, let's get the
%appropriate range indeces
%get a linear fit to the meteor's range vs time
rangeFit = polyfit((allTime(radarGoodIndeces)-allTime(radarGoodIndeces(1)))*...
    matlabTmsInSec,radarGoodRanges,1);
%use the range fit to get the range indeces where the meteor exists for the
%times we care about
radarRanges = polyval(rangeFit,(allTime(radarTimeIndeces)-...
    allTime(radarGoodIndeces(1)))*matlabTmsInSec);
% radarRanges = polyval(rangeFit,(wideBeamTime(radarTimeIndeces)-...
%     wideBeamTime(radarTimeIndeces(1)))*matlabTmsInSec);
%convert the ranges to indeces
radarRangeIndeces = zeros(1,numel(radarRanges));
for i = 1:numel(radarRangeIndeces)
   [~,rangeIndex] = min(abs(allRange - radarRanges(i)));
   radarRangeIndeces(i) = rangeIndex;
end
%now we have the time and range indeces, we can get the power of the meteor
%for the times we want!
radarPowerMeteor = zeros(1,numel(radarRangeIndeces));
radarRangeMeteor = zeros(1,numel(radarRangeIndeces));
for i = 1:numel(radarPowerMeteor)
    radarPowerMeteor(i) = allPower(radarRangeIndeces(i),radarTimeIndeces(i));
    radarRangeMeteor(i) = allRange(radarRangeIndeces(i));
end

radarTime = wideBeamTime(radarTimeIndeces) - extraTimeShift/matlabTmsInSec;
%do a polyfit to the optical data to get the times the meteor was at
%certian pixels
azCoeffs = polyfit(meteorTimes-meteorTimes(1),meteorAz,1);
elCoeffs = polyfit(meteorTimes-meteorTimes(1),meteorEl,1);
azFitRadar = polyval(azCoeffs,radarTime-meteorTimes(1));
elFitRadar = polyval(elCoeffs,radarTime-meteorTimes(1));
%get the array factor for each az and el of the radar
cd('C:\Users\Glenn\Documents\MATLAB\PFISR\')
arrayFactorsRadar = arrayFactorFromAzEl(azFitRadar,elFitRadar,beamType);
arrayFactorsOptical = arrayFactorFromAzEl(meteorAz,meteorEl,beamType);
arrayFactorMax = arrayFactorFromAzEl(15,74,beamType);
%make sure the array factor max is correct
if abs(arrayFactorMax) < max(abs(arrayFactorsRadar))
   arrayFactorMax = max(abs(arrayFactorsRadar))*1.15;
end
%get back to original dir
cd(returnDir)
%now we have the array factors
%get the meteor radar SNR
radarNoise = median(allPower(:));
radarSNRMeteor = 10*log10(radarPowerMeteor/radarNoise);
radarSNRMeteorNormalized = radarSNRMeteor + 10*log10(abs(arrayFactorMax).^2) - ...
    10*log10(abs(arrayFactorsRadar).^2);
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec, radarSNRMeteorNormalized,'r','linewidth',3)
hold on
plot((radarTime-radarTime(1))*matlabTmsInSec, radarSNRMeteor,'b','linewidth',3)
legend('Beam Normalized','Raw')
title('Meteor SNR')
xlabel('Time (sec)')
ylabel('SNR (dB)')
%set(gca,'xlim',[0 0.3])
%% END OF GETTING NORMALIZED RADAR POWER
%% CONVERT NORMALIZED RADAR POWER TO WATTS
conversionFactor = 0.1;
pRWatts = 10.^(radarSNRMeteorNormalized/10)*conversionFactor./...
    ((4*pi*(radarRangeMeteor.^2)).^2);
rcsMetersq = 10.^(radarSNRMeteorNormalized/10)*radarNoise*conversionFactor./...
    ((4*pi*(radarRangeMeteor.^2)).^2);
%get the approximate head plasma radius
cd('C:\Users\Glenn\Documents\MATLAB\Meteors')
meteorAlt = interp1(meteorTimes,meteorZ,radarTime,'spline','extrap');
meteorRadius = getPlasmaRadius(meteorAlt, velocitiesMean);
%load the spherical model
%load('C:\Users\Glenn\Documents\MATLAB\Meteors\450MHzWyatt10to1600MHz.mat')
load('C:\Users\Glenn\Documents\MATLAB\Meteors\450MHzWyatt.mat')
cd(returnDir)
%find where the meteor radius matches with the spherical model
%go through each point in time we see the meteor, and get the point on the
%spherical model we need
peakPlasmaFreqs = zeros(1,numel(radarRangeMeteor));
totalE = 0;
massLost = zeros(1,numel(radarRangeMeteor));
dEdt = zeros(1,numel(radarRangeMeteor));
for i = 1:numel(radarRangeMeteor)
    %get the rSigma Index
    [~,rSigmaIndex] = min(abs(meteorRadius(i) - rSigmasGlenn));
    rSigma = rSigmasGlenn(rSigmaIndex);
    %get the possible RCS (metersq) values
    possibleRCSmsqValues = powerDensities(:,rSigmaIndex)*4*pi;
    %find where the RCS is closest
    [~,RCSIndex] = min(abs(rcsMetersq(i) - possibleRCSmsqValues));
    %get the peak plasma freq
    peakPlasmaFreq = fpPeaks(RCSIndex);
    peakPlasmaW = peakPlasmaFreq*2*pi;
    peakPlasmaFreqs(i) = peakPlasmaFreq;
    %now we have the peak plasma frequency, so find the cross section
    %density
    plasmaDensityPeak = peakPlasmaW^2*epsilon0*me/(e^2);
    crossSectionDensity = peakPlasmaFreq*pi*rSigma^2;
    dEdt(i) = crossSectionDensity*velocitiesMean*1000*ipp;
    totalE = totalE + dEdt(i);
    if i == 1
        massLost(i) = dEdt(i);
    else
        massLost(i) = massLost(i-1) + dEdt(i);
    end
end
%get the beta
cd('C:\Users\Glenn\Documents\MATLAB\Meteors')
beta97 = getBetaJones1997(velocitiesMean);
beta01 = getBetaJones2001(velocitiesMean);
%get the mass
cd(returnDir)
totalMass97 = totalE*velocitiesMean*muEst/beta97*10;
totalMass01 = totalE*velocitiesMean*muEst/beta01*10;
totalMass = (totalMass97+totalMass01)/2;
%get the mass as a function of time
dMass97 = massLost*velocitiesMean*muEst/beta97*10;
dMass01 = massLost*velocitiesMean*muEst/beta01*10;
massLostDist = (dMass97+dMass01)/2;
%get the error in the mass:
massLostError = abs(dMass97-dMass01)/2+massLostDist*0.2;
%get the cumulative mass error
cumMassError = massLostError;
for i = 2:numel(massLostError)
    cumMassError(i) = sum(massLostError(1:i));
end

%% Plot the mass lost as a function of altitude
%get the altitude with interpolation
meteorAlt = interp1(meteorTimes,meteorZ,radarTime,'spline','extrap');
figure()
plot(meteorAlt,massLostDist)
title('Meteroid Mass')
xlabel('Alt (km)')
ylabel('Mass (kg)')

figure()
plot(meteorAlt,rcsMetersq)
title('RCS (m sq)')
xlabel('Alt (km)')
ylabel('RCS')
figure()
errorbar(meteorAlt,abs(massLostDist-massLostDist(end)),fliplr(massLostError))
%% Optical background plot
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
imagesc(backgroundSum)
set(gca,'clim',[3000 8000])
colormap gray
axis off

% %% Az, El image plots
% %load the data
% cd('C:\Users\Glenn\Documents\MATLAB\PFISR\Optical')
% load('ultraCalibrate2.mat')
% figure()
% fontsize = 20;
% az(az > 200) = az(az>200)-360;
% imagesc(az)
% set(gca,'ydir','normal')
% axis off
% colorbar
% set(gca,'fontsize',fontsize)
% title('Azimuth (deg)')
% figure()
% imagesc(el)
% set(gca,'ydir','normal')
% axis off
% colorbar
% set(gca,'fontsize',fontsize)
% title('Elevation (deg)')

%% Zoomed in Radar SNR Plot for meteor 13
figure()
fontsize = 24;
set(gca,'fontsize',fontsize)
rangeIndecesZoom = max(min(radarRangeIndeces)-20,1):...
    min(max(radarRangeIndeces)+20,numel(allRange));
imagesc((radarTime-radarTime(1))*matlabTmsInSec,rng2(rangeIndecesZoom)/1000,...
    10*log10(allPower(rangeIndecesZoom,radarTimeIndeces)/radarNoise))
set(gca,'ydir','normal')
if beamType == 1
    title('Wide Beam SNR (dB)')
else
    title('Narrow Beam SNR (dB)')
end
xlabel('Time (sec)')
ylabel('Range (km)')
colorbar
set(gca,'clim',[-3 20])
set(gca,'fontsize',fontsize)

%% Plot the Optical Data
%get the maximum SBR for each pixel
numOFrames = size(SBR,3);
SBRmax = zeros(512,512);
for i = 1:512
    for j = 1:512
        SBRmax(i,j) = max(SBR(i,j,:));
    end
end
figure()
fontsize = 24;
set(gca,'fontsize',fontsize)
imagesc(SBRmax)
colormap gray
title('Optical SBR (dB)')
axis off
colorbar
set(gca,'fontsize',fontsize)
set(gca,'clim',[0,3])
%% plot contour map of wide beam
% figure()
% set(gca,'fontsize',fontsize)
% title('Beam Pattern')
% [c,h] = contour(gainWide,[-30 -20 -10 -3],'linewidth',3);
% axis off
% title('asdf')
% colorbar
% set(gca,'fontsize',fontsize)
% clabel(c,h,'labelspacing',10000,'fontsize',fontsize,'Color','y')
% axis off

%% METEOR RCS PLOT
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec,10*log10(rcsMetersq),'linewidth',3)
title('Meteor RCS')
ylabel('Radar Cross Section (dBsm)')
xlabel('Time (sec)')

%% METEOR RADIUS PLOT
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec,meteorRadius,'linewidth',3)
title('Plasma Radius')
ylabel('Plasma Radius (m)')
xlabel('Time (sec)')

%% Meteor Mass with Errors
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
errorbar(meteorAlt,abs(massLostDist-massLostDist(end)),fliplr(massLostError),'linewidth',2)
title('Meteoroid Mass')
ylabel('Mass (kg)')
xlabel('Altitude (km)')
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
errorbar((radarTime-radarTime(1))*matlabTmsInSec,abs(massLostDist-massLostDist(end)),fliplr(massLostError),'linewidth',2)
hold on
plot((radarTime-radarTime(1))*matlabTmsInSec,abs(massLostDist-massLostDist(end)),'linewidth',3)
title('Meteoroid Mass')
ylabel('Mass (kg)')
xlabel('Time (sec)')
set(gca,'xlim',[0 0.4])

%% PLOT METEOR MASSES VS ALTITUDE FOR MULTIPLE METEORS
figure(99)
fontsize = 24;
set(gca,'fontsize',fontsize)
hold on
semilogy(meteorAlt,abs(massLostDist-massLostDist(end)),'linewidth',3)
set(gca,'xdir','reverse')
xlabel('Altitude (km)')
ylabel('Mass (kg)')
hold off
% mAlts = {mAlts{:},meteorAlt}
% mVels = {mVels{:},velocitiesMean}
% mMass = {mMass{:},massLostDist}

%% Plot all the masses
% vels = zeros(1,numel(mAlts));
% %get all the velocities
% for i = 1:numel(mAlts)
%     vels(i) = mVels{i};
% end
% maxVel = 70;
% minVel = 10;
% figure()
% fontsize = 24;
% set(gca,'fontsize',fontsize)
% for i = 2:numel(mAlts)
%     alt = mAlts{i};
%     vel = vels(i);
%     massLost = mMass{i};
%     mass = massLost(end)-massLost+massLost(1);
%     cd('C:\Users\Glenn\Documents\MATLAB\Meteors')
%     semilogy(alt,mass,'linewidth',3,'color',getRGBJet(vel,minVel,maxVel))
%     cd(returnDir)
%     set(gca,'xdir','reverse')
%     hold on
%     pause(0.3)
% end
% title('Observed Meteor Masses')
