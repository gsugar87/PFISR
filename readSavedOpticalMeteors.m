%% read in the saved meteor data and make appropriate plots
% Some useful infotmation:
% tUTCU is the times in matlab format, so datestr(tUTCU) gives the correct
% UTC time (this is in days)
% radar time (tms) is in unix time.  To convert to matlab time you need to
% do: unix_time/86400 + datenum(1970,1,1))

%some neded variables
narrowGainBoost = 40; %random boost to narrow beam gain to account for nulls
matlabTmsInSec = 86400;
%choose what meteor to read in
meteorNumber = 1;

%Set needed variables (directories, camera FPS, etc.)
directoryOptical = 'C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\Meteors\';
directoryRadar = 'C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\';
FPS = 53;
joulesPerCount = 1.2585/(6.5594E4);
tau = 0.02;                 %tau used in intensity mass calculations
radarFramesToSkip = 0;    %number of radar pulses to skip to match the times
%radarFramesToSkip = 277;    %number of radar pulses to skip to match the times

%load in the beam patterns
load('Radar\beamPatternNarrowUltra.mat')
load('Radar\beamPatternWideUltra.mat')

%read in the desired meteor radar data
switch(meteorNumber)
    case(1)
        filenameOptical = '0331082552ultra.mat';
        startFrame = 50;
        stopFrame = 160;
        maxSNROptical = 4;
        minSNROptical = -1;
        maxSNRRadar = [50 30];
        minSNRRadar = [-2 -2];
        filenameRadar = '0331082552radar.mat';
        load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\0331082552radarBefore.mat');
        odatBefore = odat;
        tmsBefore = tms;
        load([directoryRadar filenameRadar])
        odat = cat(3,odatBefore,odat);
        tms = cat(2,tmsBefore,tms);
        iend = iend*2;
        framesFromStart1 = 10;
        framesFromStart2 = 30;
        radialVelocity = (96.86 - 96.71)/(0.6405 - 0.5445);
        kmPerPixel = 4.5*pi/180/256*96.75;
    case(2)
        %frames 37 to 43
        %frame 37: (429,267)
        %frame 43: (512,320) part of it is outside, should be slightly
        %farther
        %frame 42: (499,312)
        filenameOptical = '0331112054ultra.mat';
        startFrame = 30;
        stopFrame = 105;
        maxSNROptical = 3;
        minSNROptical = -1;
        maxSNRRadar = [30 20];
        minSNRRadar = [-2 -2];
        filenameRadar = '0331112055radar.mat';
        load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\0331112055radarBefore.mat');
        odatBefore = odat;
        tmsBefore = tms;
        load([directoryRadar filenameRadar])
        odat = cat(3,odatBefore,odat);
        tms = cat(2,tmsBefore,tms);%+0.85;
        iend = iend*2;
        framesFromStart1 = 1;
        framesFromStart2 = 6;
        radialVelocity = (102.6 - 98.96)/(2.513 - 1.439);
        kmPerPixel = 4.5*pi/180/256*101;
        startCoordinates = [429, 267];
        endCoordinates = [499, 312];
        timeTravelled = (42-37)/FPS;
        distanceTravelled = sqrt((startCoordinates(1)-endCoordinates(1))^2 + ...
            (startCoordinates(2)-endCoordinates(2))^2)*kmPerPixel;
        tangentialVelocity = distanceTravelled/timeTravelled;
        totalVelocity = sqrt(radialVelocity^2 + tangentialVelocity^2);
    case(3)
        filenameOptical = '0331112929ultra.mat';
        startFrame = 15;
        stopFrame = 160;
        maxSNROptical = 12;
        minSNROptical = -1;
        maxSNRRadar = [40 25];
        minSNRRadar = [-2 -2];
        filenameRadar = '0331112930radar.mat';
        load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\0331112930radarBefore.mat');
        odatBefore = odat;
        tmsBefore = tms;
        load([directoryRadar filenameRadar])
        odat = cat(3,odatBefore,odat);
        tms = cat(2,tmsBefore,tms);
        iend = iend*2;
        framesFromStart1 = 10;
        framesFromStart2 = 30;
        radialVelocity = (86.96 - 75.72)/(2.513 - 1.346);
        kmPerPixel = 4.5*pi/180/256*80;
    otherwise
        filenameOptical = 0;
end

%load the optical data
load([directoryOptical filenameOptical])
%get the times between the radar and the optical data in the same synced
%units
timeRadar = tms(2,:) + 3600*24*datenum(1970,1,1);
timeOptical = tUTCU(startFrame:stopFrame)*3600*24;
tUTCUsec = tUTCU*3600*24;
dtOptical = timeOptical(2)-timeOptical(1); %in seconds

%get time indeces for the radar at each optical frame
opticalToRadarFrameIndeces = zeros([2,numel(timeOptical)]);
for i = 1:numel(timeOptical)
   [~,frameTo] = min(abs(timeOptical(i)-dtOptical-timeRadar));
   [~,frameTf] = min(abs(timeOptical(i)-timeRadar));
   opticalToRadarFrameIndeces(:,i) = [frameTo;frameTf];
end
opticalToRadarFrameIndeces = opticalToRadarFrameIndeces + radarFramesToSkip;
startRadarFrameIndex = min(min(opticalToRadarFrameIndeces));
stopRadarFrameIndex = max(max(opticalToRadarFrameIndeces));

dataStdU = std(double(dataUltra),[],3);
dataMeanU = mean(dataUltra,3);
signalUltra = zeros(size(dataUltra));
for i = 1:size(dataUltra,3)
   signalUltra(:,:,i) = double(dataUltra(:,:,i)) - dataUltraBackground;
end

%get the velocity in the optical data unless the meteor number is 2
%(already done)
if meteorNumber ~= 2
    %get the index of the maximum signal in a frame
    [maxValues1, maxIndeces1] = max(signalUltra);
    maxValues1 = squeeze(maxValues1);
    maxIndeces1 = squeeze(maxIndeces1);
    [maxValues2, maxIndeces2] = max(maxValues1);
    maxValues2 = squeeze(maxValues2);
    maxIndeces2 = squeeze(maxIndeces2);
    startCoordinates = [maxIndeces1(maxIndeces2(startFrame+framesFromStart1), ...
        startFrame+framesFromStart1),maxIndeces2(startFrame+framesFromStart1)];
    endCoordinates = [maxIndeces1(maxIndeces2(startFrame+framesFromStart2), ...
        startFrame+framesFromStart2),maxIndeces2(startFrame+framesFromStart2)];
    distanceTravelled = sqrt((startCoordinates(1)-endCoordinates(1))^2 + ...
        (startCoordinates(2)-endCoordinates(2))^2)*kmPerPixel;
    timeTravelled = (framesFromStart2-framesFromStart1)/FPS;
    tangentialVelocity = distanceTravelled/timeTravelled;
    totalVelocity = sqrt(radialVelocity^2 + tangentialVelocity^2);
end

%get the total signal
totalSignal = sum(sum(sum(signalUltra)));
%get the optical masses
totalJoules = totalSignal*joulesPerCount;
opticalMass = totalJoules*2/((totalVelocity*1000)^2*tau);

imagesc(dataStdU./dataMeanU)
colormap(gray)
colorbar()
title('Standard Deviation / Mean Intensity 3 sec')

figure()
imagesc(10*log10(dataMeanU./dataUltraBackground))
colormap(gray)
colorbar()
title('Signal to Background Ratio')

figure()
imagesc(10*log10(dataMeanU./median(median(dataUltraBackground))))
colormap(gray)
colorbar()
title('Signal to Noise Ratio')

figure()
%radar data
ibms = cellstr(['Narrow';'Wide  ']);
for ibm = 1:2
    x = squeeze(odat(ibm,:,ibeg:iend));
    
    subplot(2,1,ibm);
    imagesc(tms(ibm,ibeg:iend)-tm(1,irec),rng2/1e3,10*log10(x./median(median(x))));
    set(gca,'ydir','normal');
    %ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
    ylabel(['Range (km);' ibms(ibm) ' Beam'])
    colorbar;
end

%JOINT MOVIE
%makeJointMovie

%Get intensity vs time for optical data
%Track the meteor, get its coordinates as a function of time
%go to the middle of the image and find where the meteor is
numOpticalFrames = size(dataUltra,3);
maxXs = zeros(1,numOpticalFrames);
maxYs = zeros(1,numOpticalFrames);
intensityMaxsDB = zeros(1,numOpticalFrames);
intensityMaxsNarrowNormDB = zeros(1,numOpticalFrames);
intensityMaxsWideNormDB = zeros(1,numOpticalFrames);
intensityTimes = linspace(0,(numOpticalFrames-1)/FPS,numOpticalFrames);
for i = 1:numOpticalFrames
    frameSBR = 10*log10(double(dataUltra(:,:,i))./dataUltraBackground);
    [frameSBRMax,frameSBRMaxY] = max(max(frameSBR));
    frameSBRMaxX = mod(find(frameSBR == frameSBRMax),size(dataUltra,1));
    if frameSBRMaxX == 0
        frameSBRMaxX = 1;
    end
    maxXs(i) = frameSBRMaxX;
    maxYs(i) = frameSBRMaxY;
    intensityMaxsDB(i) = 10*log10(double(dataUltra(frameSBRMaxX,frameSBRMaxY,i)));
    intensityMaxsNarrowNormDB(i) = intensityMaxsDB(i)+gainNarrow(frameSBRMaxX,frameSBRMaxY);
    intensityMaxsWideNormDB(i) = intensityMaxsDB(i)+gainWide(frameSBRMaxX,frameSBRMaxY);
end
%get the frames we actually see the meteor well
possibleGoodFrames = find((abs(maxXs(2:end)-maxXs(1:end-1)) < 30) & (abs(maxYs(2:end)-maxYs(1:end-1)) < 30));
possibleGoodFramesCheck = (possibleGoodFrames(2:end)-possibleGoodFrames(1:end-1))==1;
goodFrames = possibleGoodFrames([possibleGoodFramesCheck possibleGoodFramesCheck(end)]);
%adjust the normalized intensities so that convolution will work
iMaxsNarNormDBZeroed = intensityMaxsNarrowNormDB - min(intensityMaxsNarrowNormDB(goodFrames));
iMaxsWidNormDBZeroed = intensityMaxsWideNormDB - min(intensityMaxsWideNormDB(goodFrames));
%plot the optical intensities normalized to the beam pattern
figure()
plot(tUTCUsec(goodFrames)-tUTCUsec(goodFrames(1)),intensityMaxsNarrowNormDB(goodFrames))
title('Optical Max Intensity (dB) with Normalized Values')
xlabel('Time (sec)')
ylabel('Maximum Intensity (dB)')
hold on
plot(tUTCUsec(goodFrames)-tUTCUsec(goodFrames(1)),intensityMaxsWideNormDB(goodFrames),'r')
plot(tUTCUsec(goodFrames)-tUTCUsec(goodFrames(1)),intensityMaxsDB(goodFrames),'g')
legend('Narrow Beam','Wide Beam','No Normalization')
hold off

%get the power profile for the meteor in the radar data
powerNarrow = squeeze(odat(1,:,ibeg:iend));
powerNarrowDB = 10*log10(powerNarrow./median(median(powerNarrow)));
powerNarrowDBTimes = tms(1,ibeg:iend);
powerNarrowDBMaxs = max(powerNarrowDB);
powerWide = squeeze(odat(2,:,ibeg:iend));
powerWideDB = 10*log10(powerWide./median(median(powerWide)));
powerWideDBTimes = tms(2,ibeg:iend);
powerWideDBMaxs = max(powerWideDB);

figure()
plot(powerNarrowDBTimes-powerNarrowDBTimes(1),powerNarrowDBMaxs)
title('Maximum Power (dB)')
xlabel('Time (sec)')
ylabel('Max Power (dB)')
hold on
plot(powerWideDBTimes-powerWideDBTimes(1),powerWideDBMaxs,'r')
legend('Narrow Beam','Wide Beam')
hold off

%overplot the radar signal with the optical intensity
%convert the radar times (unix) to matlab times
% note that the narrow beam times are the first dimension of tms
tmsMatlab = tms/matlabTmsInSec + datenum(1970,1,1);
tmsMatNar = tmsMatlab(1,:);
tmsMatWid = tmsMatlab(2,:);
dtMatlab = mean(tmsMatlab(1,2:end)-tmsMatlab(1,1:end-1));
figure()
plot(tmsMatNar(ibeg:iend),powerNarrowDBMaxs)
datetick('x','HH:MM:SS')
xlabel(['UTC Time On ', datestr(tmsMatlab(1),'mm/dd/yyyy')])
ylabel('Power or Intensity (dB)')
title('Radar and Optical Signals Before Time Shift')
hold on
plot(tmsMatWid(ibeg:iend),powerWideDBMaxs,'r')
plot(tUTCU(goodFrames),intensityMaxsNarrowNormDB(goodFrames),'g')
plot(tUTCU(goodFrames),intensityMaxsWideNormDB(goodFrames),'c')
plot(tUTCU(goodFrames),intensityMaxsDB(goodFrames),'k')
legend('Radar Narrow','Radar Wide','Camera Narrow','Camera Wide','Camera')

%do a convolution to try and find the time shift
%first we need to get the optical data on the same time as the radar data
%for the convolution, so we will use interpolation
%iMaxNarrowInterpNans = interp1(tUTCU(goodFrames),intensityMaxsNarrowNormDB(goodFrames),tmsMatNar);
%iMaxWideInterpNans = interp1(tUTCU(goodFrames),intensityMaxsWideNormDB(goodFrames),tmsMatWid);
iMaxNarrowInterpNans = interp1(tUTCU(goodFrames),iMaxsNarNormDBZeroed(goodFrames),tmsMatWid);
iMaxWideInterpNans = interp1(tUTCU(goodFrames),iMaxsWidNormDBZeroed(goodFrames),tmsMatWid);
iMaxNarrowInterp = iMaxNarrowInterpNans(isfinite(iMaxNarrowInterpNans));
iMaxWideInterp = iMaxWideInterpNans(isfinite(iMaxWideInterpNans));
iMaxNarrowTimes = tmsMatNar(isfinite(iMaxNarrowInterpNans));
iMaxWideTimes = tmsMatWid(isfinite(iMaxWideInterpNans));

iMaxNarrowInterp = iMaxNarrowInterp - min(iMaxNarrowInterp);
iMaxWideInterp = iMaxWideInterp - min(iMaxWideInterp);

%plot the interpolated data to make sure it worked
figure()
plot((tmsMatNar-tmsMatNar(1))*24*3600,powerNarrowDBMaxs)
hold on
plot((tmsMatWid-tmsMatWid(1))*24*3600,powerWideDBMaxs,'r')
plot((iMaxNarrowTimes-tmsMatNar(1))*24*3600,iMaxNarrowInterp,'g')
plot((iMaxWideTimes-tmsMatWid(1))*24*3600,iMaxWideInterp,'k')
title('Time Interpolated Optical Signals Before Shift')
xlabel('Time (Sec)')
ylabel('Intensity or Power (dB)')
legend('Radar Narrow','Radar Wide','Camera Narrow','Camera Wide')

%do a convolution to find the time shift between the signals
convWide = conv(fliplr(powerWideDBMaxs),iMaxWideInterp);
convNarrow = conv(fliplr(powerNarrowDBMaxs),iMaxNarrowInterp);

[~,timeShiftIndexWide] = max(abs(convWide));
[~,timeShiftIndexNarrow] = max(abs(convNarrow));
%now timeShiftIndex is the number of indeces to shift the optical signal
%forward (or the radar signal backward) to get the best alignment
timeShiftIndexWide = length(powerWideDBMaxs)-timeShiftIndexWide;
timeShiftIndexNarrow = length(powerNarrowDBMaxs)-timeShiftIndexNarrow;
%calculate the time shift
%get the time shift from the convoultion
timeShiftConvNarrow = timeShiftIndexNarrow*dtMatlab - (iMaxNarrowTimes(1)-tmsMatNar(1));
timeShiftConvWide = timeShiftIndexWide*dtMatlab - (iMaxWideTimes(1)-tmsMatWid(1));
%timeShiftConvNarrow = timeShiftConvWide;

%plot the data with the time shift of the radar
figure()
plot(tmsMatlab(1,ibeg:iend)-timeShiftConvNarrow,powerNarrowDBMaxs)
%datetick('x','HH:MM:SS')
xlabel(['UTC Time On ', datestr(tmsMatlab(1),'mm/dd/yyyy')])
ylabel('Power or Intensity (dB)')
title('Radar and Optical Signals After Time Shift')
hold on
plot(tmsMatlab(2,ibeg:iend)-timeShiftConvWide,powerWideDBMaxs,'r')
plot(tUTCU(goodFrames),iMaxsNarNormDBZeroed(goodFrames),'g')
plot(tUTCU(goodFrames),iMaxsWidNormDBZeroed(goodFrames),'c')
plot(tUTCU(goodFrames),intensityMaxsDB(goodFrames),'k')
legend('Radar Narrow','Radar Wide','Camera Narrow','Camera Wide','Camera')
hold off

middleFrameNum = round(numOpticalFrames/2);
middleFrameSBR = 10*log10(double(dataUltra(:,:,middleFrameNum))./dataUltraBackground);
[middleFrameSBRMax,middleFrameSBRMaxY] = max(max(middleFrameSBR));
middleFrameSBRMaxX = mod(find(middleFrameSBR == middleFrameSBRMax),size(dataUltra,1));
%go backwards in time to get the 
