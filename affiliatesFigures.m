%% PLOTS FOR AFFILIATES POSTER %%
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

extraTimeShift = 0.05;
%% GET THE NORMALIZED RADAR POWER %%
beamType = 1; %assume wide beam
%first we need to expand the "good indeces" of the radar for the meteor
% meteor 9 (startTime = 3.93 endTime = 4.23)
%1: 4.6, 4.8
%2: 3.32, 3.6
%3: 1.7, 2
%5: 2.35, 2.95
%6: 1.3, 1.9
%7: 1.9, 2.05
startTime = 1.72;
endTime = 2.02;
[~,radarStartIndexTime] = min(abs((wideBeamTime-wideBeamTime(1))*matlabTmsInSec - startTime));
[~,radarEndIndexTime] = min(abs((wideBeamTime-wideBeamTime(1))*matlabTmsInSec - endTime));
%get the time indeces that the meteor exists
radarTimeIndeces = linspace(radarStartIndexTime,...
    radarEndIndexTime, radarEndIndexTime-radarStartIndexTime + 1);

%now that we have the time indeces we want to plot, let's get the
%appropriate range indeces
%get a linear fit to the meteor's range vs time
rangeFit = polyfit((wideBeamTime(radarGoodIndeces)-wideBeamTime(radarGoodIndeces(1)))*...
    matlabTmsInSec,radarGoodRanges,1);
%use the range fit to get the range indeces where the meteor exists for the
%times we care about
radarRanges = polyval(rangeFit,(wideBeamTime(radarTimeIndeces)-...
    wideBeamTime(radarTimeIndeces(1)))*matlabTmsInSec);
%convert the ranges to indeces
radarRangeIndeces = zeros(1,numel(radarRanges));
for i = 1:numel(radarRangeIndeces)
   [~,rangeIndex] = min(abs(wideBeamRange - radarRanges(i)));
   radarRangeIndeces(i) = rangeIndex;
end
%now we have the time and range indeces, we can get the power of the meteor
%for the times we want!
radarPowerMeteor = zeros(1,numel(radarRangeIndeces));
radarRangeMeteor = zeros(1,numel(radarRangeIndeces));
for i = 1:numel(radarPowerMeteor)
    radarPowerMeteor(i) = wideBeamPower(radarRangeIndeces(i),radarTimeIndeces(i));
    radarRangeMeteor(i) = wideBeamRange(radarRangeIndeces(i));
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
%get back to original dir
cd(returnDir)
%now we have the array factors
%get the meteor radar SNR
radarNoise = median(wideBeamPower(:));
radarSNRMeteor = 10*log10(radarPowerMeteor/radarNoise);
radarSNRMeteorNormalized = radarSNRMeteor + 10*log10(abs(arrayFactorMax).^2) - ...
    10*log10(abs(arrayFactorsRadar).^2);
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec, radarSNRMeteorNormalized,'linewidth',3)
hold on
plot((radarTime-radarTime(1))*matlabTmsInSec, radarSNRMeteor,'r','linewidth',3)
legend('Beam Normalized','Raw')
title('Meteor SNR')
xlabel('Time (sec)')
ylabel('SNR (dB)')
set(gca,'xlim',[0 0.3])
%% END OF GETTING NORMALIZED RADAR POWER
%% CONVERT NORMALIZED RADAR POWER TO WATTS
conversionFactor = 0.1;
pRWatts = 10.^(radarSNRMeteorNormalized/10)*conversionFactor./...
    ((4*pi*(radarRangeMeteor.^2)).^2);
rcsMetersq = 10.^(radarSNRMeteorNormalized/10)*radarNoise*conversionFactor./...
    ((4*pi*(radarRangeMeteor.^2)).^2);
%get the approximate head plasma radius
cd('C:\Users\Glenn\Documents\MATLAB\Meteors')
meteorRadius = getPlasmaRadius(radarRangeMeteor, velocitiesMean);
%load the spherical model
load('C:\Users\Glenn\Documents\MATLAB\Meteors\450MHzWyatt.mat')
cd(returnDir)
%find where the meteor radius matches with the spherical model
%go through each point in time we see the meteor, and get the point on the
%spherical model we need
peakPlasmaFreqs = zeros(1,numel(radarRangeMeteor));
totalE = 0;
dmass = zeros(1,numel(radarRangeMeteor));
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
    totalE = totalE + crossSectionDensity*velocitiesMean*1000*ipp;
    dmass(i) = crossSectionDensity*velocitiesMean*1000*ipp;
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

%% Plot the radius vs time
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec,meteorRadius*0.023,'linewidth',3)
set(gca,'xlim',[0 0.3])
xlabel('Time (sec)')
ylabel('Plasma Radius (m)')
title('Plasma Radius')

%% Plot the RCS vs time
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot((radarTime-radarTime(1))*matlabTmsInSec, 10*log10(rcsMetersq),'linewidth',3)
set(gca,'xlim',[0 0.3])
title('Meteor RCS')
xlabel('Time (sec)')
ylabel('Radar Cross Section (dBsm)')


%plot the radar data
figure(101)
fontsize = 20;
movegui('northwest')
%subplot(2,1,1)
set(gca,'fontsize',fontsize)
imagesc((wideBeamTime-wideBeamTime(1))*matlabTmsInSec,wideBeamRange,10*log10(wideBeamPower/median(wideBeamPower(:))))
%imagesc((wideBeamTime-wideBeamTime(1))*matlabTmsInSec-1.7,wideBeamRange,10*log10(wideBeamPower/median(wideBeamPower(:))))
xlabel('Time (sec)')
ylabel('Range (km)')
title('Radar SNR (dB)')
set(gca,'ydir','normal')
set(gca,'xlim',[0 0.3])
set(gca,'ylim',[75 115])
colorbar
set(gca,'fontsize',fontsize)
figure(100)
%subplot(2,1,2)
set(gca,'fontsize',fontsize)
imagesc((narrowBeamTime-wideBeamTime(1))*matlabTmsInSec-1.7,narrowBeamRange,10*log10(narrowBeamPower/median(narrowBeamPower(:))))
set(gca,'ydir','normal')
set(gca,'xlim',[0 1])
set(gca,'ylim',[70 110])
xlabel('Time (sec)')
ylabel('Range (km)')
title('Narrow Beam SNR (dB)')
colorbar

%%OPTICAL IMAGE %%
figure()
set(gca,'fontsize',fontsize)
imagesc(max(SBR,[],3))
title('Optical SBR (dB)')
axis off
colormap gray
colorbar
set(gca,'fontsize',fontsize)

%%OPTICAL IMAGE %%
figure()
set(gca,'fontsize',fontsize)
imagesc(dataUltraBackground)
axis off
colormap gray
colorbar
set(gca,'fontsize',fontsize)

%%Beam Pattern%%
figure()
%fontsize = 18;
load('beamPatternWideUltra.mat')
[c,h] = contour(gainWide,[-30 -20 -10 -3],'linewidth',3);
title('Beam Pattern','fontsize',20)
colormap jet
clabel(c,h,'labelspacing',10000,'fontsize',fontsize,'Color','y')
colorbar
set(gca,'fontsize',fontsize)
axis off

%%Spherical Modeling Plot
load('C:\Users\Glenn\Documents\MATLAB\Meteors\450MHzWyatt.mat')
fontsize = 20;
rcs = 10*log10(powerDensities*4*pi);
imagesc(rSigmasGlenn,fpPeaks/1e6,rcs)
set(gca,'fontsize',fontsize)
%set(gca,'xlim',[0,5])
set(gca,'ydir','normal')
set(gca,'clim',[-60,5])
colorbar
set(gca,'fontsize',fontsize)
title('RCS (dBsm)')
xlabel('Head Plasma Radius (m)')
ylabel('Peak Plasma Freq (MHz)')

%% ALTITUDE VELOCITY MASS PLOT
masses = [0.00000000029985;
          0.0000000027655;
          0.000000017171;
          0.0000000084544;
          0.0000000054044;
          0.00000000055026;
          0.000000035867];
vels = [28.2757;
        60.3688;
        28.0872;
        58.6723;
        22.6498;
        22.8142;
        18.7481];
maxAlts = [94.5661;
98.4309;
99.0238;
95.2066;
102.1429;
83.8601;
99.86];
minAlts = [86.2255;
95.7711;
94.759;
91.7478;
82.2832;
79.9567;
89.0679];
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
for i = 1:numel(minAlts)
    plot([vels(i),vels(i)],[maxAlts(i),minAlts(i)],'color',...
        getRGBJet(10*log10(masses(i)),10*log10(min(masses)),...
        10*log10(max(masses))),'linewidth',4)
    hold on
end
title('Observed Meteors')
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')
set(gca,'ylim',[75 105])

%% Plot gaussian distribution
rDist = linspace(0,3,1000);
wp = 500*exp(-(rDist/1).^2);
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
plot(rDist,wp,'linewidth',3)
xlabel('Radius (m)')
ylabel('Plasma Freq (MHz)')
title('\sigma_r = 1 m, \omega_0 = 500 MHz')

%% Plot plasma radius
v = ones(1,1000)*25;
alts = linspace(80,110,1000);
cd('C:\Users\Glenn\Documents\MATLAB\Meteors')
rs = getPlasmaRadius(alts, v);
cd(returnDir)
figure()
fontsize = 20;
set(gca,'fontsize',fontsize)
semilogx(rs,alts,'linewidth',3)
xlabel('\sigma_r (m)')
ylabel('Altitude (km)')
title('Plasma Radius for v = 25 km/s')
