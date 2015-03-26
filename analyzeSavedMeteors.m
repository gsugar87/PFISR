%% Updated read in the saved meteor data and do appropriate analysis
% Some useful infotmation:
% tUTCU is the times in matlab format, so datestr(tUTCU) gives the correct
% UTC time (this is in days)
% radar time (tms) is in unix time.  To convert to matlab time you need to
% do: unix_time/86400 + datenum(1970,1,1))
clear all
close all

%choose what meteor to read in
for meteorNumber = [1 2 3 4 5 6 7 9];
%meteorNumber = 4;

%set show movies or not
show_SBR = 1; %show a movie of the signal to background ratio
show_signal = 1; %%show a movie of the signal (intensity - background)
plotPauseTime = 0.01; %pause time between SBR strip plots

%set some needed constants
opticalToRadarTimeShiftSec = 1;
matlabTmsInSec = 86400;
directoryOptical = 'C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\Meteors\';
directoryRadar = 'C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\';
FPS = 53;
joulesPerCount = 1.2585/(6.5594E4); %joulse per pixel count in optical data
tau = 0.02;                 %tau used in intensity mass calculations
radarFramesToSkip = 0;    %number of radar pulses to skip to match the times
meteorPixelThickness = 2; %number of pixels to look up and down that approximates the meteor thickness in optical data in pixels
%optical masking and meteor determination variables
allowedBlankPoints = 1; %this is the allowed blank points for checking the end of a meteor in the optical data
SBRMaskthresh = 1.4; %1.2;
SBRsmallerThresh = 0.6;
SBREndPointThresh = 0.6; %0.6
pixelVelThreshFactor = 0.2; %factor to throw away possible meteor points in the optical frame
pRangeFit = 1; %degree of polynomial to fit the radar range for optical range extrapolation
%rcs variables
powerTransmit = 2e6; %2 megaWatts (check this!)
radarGain = 43; %43dbi gain
freq = 445;
c = 2.998e8;
wavelength = c/freq;

%load in the beam patterns
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\beamPatternNarrowUltra.mat')
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\beamPatternWideUltra.mat')

%create the filenames
opticalFilenames = {'0330105633ultra.mat',
    '0330105754ultra.mat',
    '0330110015ultra.mat',
    '0330111636ultra.mat',
    '0330112439ultra.mat',
    '0330114821ultra.mat',
    '0330121306ultra.mat',
    '0330122329ultra.mat',
    '0330123903ultra.mat'};
radarFilenames = {'0331105634.mat',
    '0331105754.mat',
    '033111015.mat',
    '0331111637.mat',
    '0331112440.mat',
    '0331114822.mat',
    '033112137.mat',
    '0331122329.mat',
    '033112394.mat'};
totalFramesAllMeteors = [13, 26, 14, 34, 55, 25, 41, 19, 17, 25];
startFramesAllMeteors = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20];
goodHoughIndeces = {[2 3 4 5];
    [2 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ]};

%read in the desired meteor
filenameOptical = opticalFilenames{meteorNumber};
filenameRadar = radarFilenames{meteorNumber};
numOpticalFrames = totalFramesAllMeteors(meteorNumber);
startFrame = startFramesAllMeteors(meteorNumber);

%load the radar and optical data
load([directoryRadar filenameRadar])
load([directoryOptical filenameOptical])

%% Optical Analysis
%convert data type to double and extract the data we care about
dataUltra = double(dataUltra(:,:,startFrame:startFrame+numOpticalFrames));
dataUltraBackground = double(dataUltraBackground);
tUTCU = tUTCU(startFrame:startFrame+numOpticalFrames);

%get the SBR (signal to background noise)
SBR = zeros(size(dataUltra));
maxSBR = zeros(size(dataUltra,3));
for i = 1:numOpticalFrames+1
    SBR(:,:,i) = 10*log10(dataUltra(:,:,i)./dataUltraBackground);
    maxSBR(i) = max(max(SBR(:,:,i)));
end

%get other important data from the data
dataStd = std(dataUltra,[],3);
dataMean = mean(dataUltra,3);
dataSum = sum(dataUltra,3);
backgroundSum = dataUltraBackground*size(dataUltra,3);

%make the times between the radar and the optical data in the same synced
%units
timeRadar = tms(2,:) + 3600*24*datenum(1970,1,1);
timeOptical = tUTCU*3600*24;
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

%get the signal (intensity - background) for the ultra data
dataStdU = std(double(dataUltra),[],3);
dataMeanU = mean(dataUltra,3);
signalUltra = zeros(size(dataUltra));
for i = 1:size(dataUltra,3)
    signalUltra(:,:,i) = double(dataUltra(:,:,i)) - dataUltraBackground;
end

%% Find the streaks in the image using a Hough Transform
% NEW METHOD ON 3/9/2015
%Step 1) mask each SBR frame, then add them together to get a masked total
%frame
SBRmaskedAll = zeros(size(SBR));
for i = 1:size(SBR,3)
    %mask the images
    SBRmaskBuffer = zeros(size(SBR,1),size(SBR,2));
    SBRmaskBuffer(SBR(:,:,i) >= SBRMaskthresh) = 1;
    SBRmaskedAll(:,:,i) = SBRmaskBuffer;
end
SBRmasked = sum(SBRmaskedAll,3);
SBRmasked = SBRmasked > 0;
[H,theta,rho] = hough(SBRmasked);
%find rho and theta of maximum
[HMaxCol,maxCol] = max(H);
[~,maxRow] = max(HMaxCol);
rhoMax = rho(maxCol(maxRow));
thetaMax = theta(maxRow);
%see if we need to do a circle mask
if abs(thetaMax) == 45
    %we need to do a circle mask
    %circle mask
    [xCircle,yCircle] = meshgrid(-255:256,-255:256);
    circleMask = ((xCircle.^2+yCircle.^2)<=256^2);
    [H,theta,rho] = hough(SBRmasked.*circleMask);
    %find rho and theta of maximum
    [HMaxCol,maxCol] = max(H);
    [~,maxRow] = max(HMaxCol);
    rhoMax = rho(maxCol(maxRow));
    thetaMax = theta(maxRow);
end
wideBeamPower = cosd(thetaMax)*rhoMax;
y = sind(thetaMax)*rhoMax;
slope = -1/tand(thetaMax);
intercept = y - slope*wideBeamPower;
yCoords = zeros(1,size(SBR,3)).*nan;
xCoords = zeros(1,size(SBR,3)).*nan;
%(TESTING) Plot the line on the net SBRmasked image
figure(1)
movegui('northwest')
imagesc(SBRmasked)
hold on
plot(linspace(1,512,512),linspace(1,512,512)*slope+intercept,'g')
title('Masked SBR of All Frames with Hough Line')
hold off
%(END TESTING)
%Now that we have the line the meteor is on, go through each frame and find
%the beginning and end of the meteor
%Do this by going through each x index (if thetaMax > 45) or y index 
%(if thetaMax > 45) and find the
%biggest batch of positive masked values and get the beginning and end
%points for it
SBRstrips = zeros(size(SBR,3),size(SBR,1));
SBRstripsX = zeros(size(SBR,3),size(SBR,1));
SBRstripsY = zeros(size(SBR,3),size(SBR,1));
for i = 1:size(SBR,3)
   %go through each image
   for j = 1:size(SBR,1)
       %go through each x and y index
       yUp = round(intercept+j*slope+meteorPixelThickness);
       yDown = round(intercept+j*slope-meteorPixelThickness);
       xUp = round((j-intercept)/slope+meteorPixelThickness);
       xDown = round((j-intercept)/slope-meteorPixelThickness);
       if abs(thetaMax) >= 45
           %we will use the x coodinate to find the appropriate y coords
           if yUp > 0 && yDown < size(SBR,1)
               SBRstrips(i,j) = mean(SBR(max([yDown 1]):min([yUp size(SBR,1)]),j,i));
           end
       else
           if xUp > 0 && xDown < size(SBR,1)
               SBRstrips(i,j) = mean(SBR(j,max([xDown 1]):min([xUp size(SBR,1)]),i));
           end
       end
   end
end
% end of going through each index to get the SBR strip
%see if the meteor is travelling forwards or backwards
[~,SBRstripsMaxIndeces] = max(SBRstrips');
if median(SBRstripsMaxIndeces(2:end)-SBRstripsMaxIndeces(1:end-1)) >= 0
    %the meteor is travelling "forwards"
    travelDir = 1;
else
    %meteor is travelling "backwards"
    travelDir = -1;
end
%go through each image again
for i = 1:size(SBR,3)
   %look at the SBR strips for the i'th image and find the longest streak
   %create an array to hold possible meteors [length,start,end, SBRmean]
    possibleMeteors = [[],[],[],[]];
    j = 1;
    if travelDir == -1
        SBRstripLook = fliplr(SBRstrips(i,:));
    else
        SBRstripLook = SBRstrips(i,:);
    end
    while j < size(SBR,1)
        %see if the index is a potential meteor
        if SBRstripLook(j) > 0
            %there is a point over the SBR threshold
            %see if there is another point next to it that is over the SBR
            %threshold
            validPoints = find(SBRstripLook(j+1:min([j+1+allowedBlankPoints size(SBR,1)])) > SBREndPointThresh);
            if numel(validPoints) > 0
                %this is a potential meteor
                potentialMeteorStartIndex = j;
                potentialMeteorEndIndex = j+validPoints(end);
                %find the end point of the potential meteor
                while potentialMeteorEndIndex <= size(SBR,1) && numel(validPoints) > 0
                    validPoints = find(SBRstrips(i,potentialMeteorEndIndex+1:min([potentialMeteorEndIndex+1+allowedBlankPoints size(SBR,1)])) > 0);
                    if numel(validPoints) > 0
                        potentialMeteorEndIndex = validPoints(end)+potentialMeteorEndIndex;
                    end
                end
                potentialMeteorMeanSBR = sum(SBRstripLook(potentialMeteorStartIndex:potentialMeteorEndIndex));
                %save the potential meteor
                if travelDir == 1
                    possibleMeteors = [possibleMeteors; potentialMeteorEndIndex-potentialMeteorStartIndex+1 potentialMeteorStartIndex potentialMeteorEndIndex potentialMeteorMeanSBR];
                else
                    possibleMeteors = [possibleMeteors; potentialMeteorEndIndex-potentialMeteorStartIndex+1 size(SBR,1)-potentialMeteorStartIndex+1 size(SBR,1)-potentialMeteorEndIndex+1 potentialMeteorMeanSBR];
                end
                %increment j
                j = potentialMeteorEndIndex+1;
            else
                j = j+1;
            end
        else
            %increase j by 1
            j = j+1;
        end
    end
    %get the max and min coordinates for the possible meteor with the
    %largest mean SBR
    if numel(possibleMeteors) > 0
        [~,possibleMeteorIndex] = max(possibleMeteors(:,4));
        if abs(thetaMax) >= 45
            %we use the x coords and find the y
            xCoords(i) = possibleMeteors(possibleMeteorIndex,3);
            yCoords(i) = round(intercept + xCoords(i)*slope);
        else 
            yCoords(i) = possibleMeteors(possibleMeteorIndex,3);
            xCoords(i) = round((yCoords(i)-intercept)/slope);
        end
    end
end



%(TESTING) do a time lapse plot of the coordinates
for i = 1:size(SBR,3)
    figure(2)
    movegui('north')
    plot(SBRstrips(i,:))
    hold on
    if abs(thetaMax) >= 45
        plot([xCoords(i) xCoords(i)],[min(SBRstrips(i,:)) max(SBRstrips(i,:))],'r')
    else
        plot([yCoords(i) yCoords(i)],[min(SBRstrips(i,:)) max(SBRstrips(i,:))],'r')
    end
    title(['SBR along Meteor Line for Frame ' num2str(i)])
    xlabel('X pixel')
    ylabel('SBR')
    hold off
    figure(3)
    movegui('northeast')
    hold on
    set(gca,'xlim',[0 512])
    set(gca,'ylim',[0 512])
    plot(xCoords(i),yCoords(i),'x')
    text(xCoords(i),yCoords(i), num2str(i), 'VerticalAlignment','bottom', ...
                                            'HorizontalAlignment','right')
    title('The Pixel Position of the Meteor at Each Frame')
    hold off
    pause(plotPauseTime)
end
%(END TESTING)
%we have the potential x and y coordinates of the meteor, but there could
%be outliers (points that are not the meteor)
%find all the valid frames (not nans)
validFrames = find(isfinite(xCoords) == 1);
frameNum = linspace(1,numel(xCoords),numel(xCoords));
%get the pixels travelled at each frame
%look at the pixel distance from the middle valid frame
middleValidFrameCoord = [xCoords(validFrames(ceil(end/2))) yCoords(validFrames(ceil(end/2)))];
middleValidFrameNum = validFrames(ceil(end/2));
pixelVelocities = sqrt((xCoords - middleValidFrameCoord(1)).^2 + (yCoords - middleValidFrameCoord(2)).^2)./abs(frameNum - middleValidFrameNum);
pixelVelocities(middleValidFrameNum) = nanmedian(pixelVelocities);
%find where pixel velocities are within a set amount of the median velocity
goodFrames = find(abs(pixelVelocities - nanmedian(pixelVelocities)) < nanmedian(pixelVelocities)*pixelVelThreshFactor);
%goodFrames = abs(pixelsTravelled - nanmedian(pixelsTravelled)) < 2*nanstd(pixelsTravelled);
%now that we have good frames and the coordinates, we need to find the
%azimuth and elevation of the meteor in each frame
%read in the ultra calibration
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\ultraCalibrate2.mat')
%read in beam patterns
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\beamPatternWideUltra.mat')
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\beamPatternNarrowUltra.mat')
%first adjust for possible looping
az(az > 150) = az(az > 150) - 360;
%create arrays to hold meteor information
meteorAz = zeros(1,numel(goodFrames));
meteorEl = zeros(1,numel(goodFrames));
meteorFrame = zeros(1,numel(goodFrames));
meteorTimes = zeros(1,numel(goodFrames));
meteorWideGain = zeros(1,numel(goodFrames));
meteorNarrowGain = zeros(1,numel(goodFrames));
for i = 1:numel(goodFrames)
    meteorAz(i) = az(yCoords(goodFrames(i)),xCoords(goodFrames(i)));
    meteorEl(i) = el(yCoords(goodFrames(i)),xCoords(goodFrames(i)));
    meteorFrame(i) = goodFrames(i); %the frame of the optical data
    meteorTimes(i) = tUTCU(goodFrames(i)); %the optical data time for the frame
    meteorWideGain(i) = gainWide(yCoords(goodFrames(i)),xCoords(goodFrames(i)));
    meteorNarrowGain(i) = gainNarrow(yCoords(goodFrames(i)),xCoords(goodFrames(i)));
end
%Now we have the azimuth, elevation, and time for the meteor at valid
%frames
%Now we need to use the radar data to get the range of the meteor at the
%time of the valid frames
%get the wide beam information
wideBeamPower = squeeze(odat(2,:,ibeg:iend));
wideBeamTime = tms(2,ibeg:iend)./matlabTmsInSec + datenum(1970,1,1) - opticalToRadarTimeShiftSec/(24*60*60);
narrowBeamPower = squeeze(odat(1,:,ibeg:iend));
narrowBeamTime = tms(1,ibeg:iend)./matlabTmsInSec + datenum(1970,1,1) - opticalToRadarTimeShiftSec/(24*60*60);
wideBeamRange = rng2/1e3;
narrowBeamRange = rng2/1e3;
%find out where the meteor is in the wide beam
%look at the net power at all times
wideBeamNetPower = sum(wideBeamPower);
wideBeamValidPoints = wideBeamNetPower > 1.1*nanmedian(wideBeamNetPower);
narrowBeamNetPower = sum(narrowBeamPower);
narrowBeamValidPoints = narrowBeamNetPower > 1.1*nanmedian(narrowBeamNetPower);
%find the longest streak of valid points (this is the meteor)
possibleMeteorStart = 1;
possibleMeteorEnd = 1;
possibleMeteorLength = 1;
possibleMeteorStartBest = 1;
possibleMeteorEndBest = 1;
possibleMeteorLengthBest = 1;
for i = 2:numel(wideBeamNetPower)
    if wideBeamValidPoints(i) == 1
        %see if the previous point was also valid
        if wideBeamValidPoints(i-1) == 1
            %increase the possible length
            possibleMeteorLength = possibleMeteorLength +1;
        else
            %get a new starting point
            possibleMeteorStart = i;
        end
    else
        %see if the previous point was valid
        if wideBeamValidPoints(i-1) == 1
            %update the end point
            possibleMeteorEnd = i-1;
            possibleMeteorLength = possibleMeteorEnd-possibleMeteorStart;
            %see if the length of this meteor is better than the best
            if possibleMeteorLength > possibleMeteorLengthBest
                %update the best possible meteor
                possibleMeteorStartBest = possibleMeteorStart;
                possibleMeteorEndBest = possibleMeteorEnd;
                possibleMeteorLengthBest = possibleMeteorLength;
            end
        end
    end
end
%do the same for the narrow beam
possibleMeteorStartNarrow = 1;
possibleMeteorEndNarrow = 1;
possibleMeteorLengthNarrow = 1;
possibleMeteorStartBestNarrow = 1;
possibleMeteorEndBestNarrow = 1;
possibleMeteorLengthBestNarrow = 1;
for i = 2:numel(narrowBeamNetPower)
    if narrowBeamValidPoints(i) == 1
        %see if the previous point was also valid
        if narrowBeamValidPoints(i-1) == 1
            %increase the possible length
            possibleMeteorLengthNarrow = possibleMeteorLengthNarrow +1;
        else
            %get a new starting point
            possibleMeteorStartNarrow = i;
        end
    else
        %see if the previous point was valid
        if narrowBeamValidPoints(i-1) == 1
            %update the end point
            possibleMeteorEndNarrow = i-1;
            possibleMeteorLengthNarrow = possibleMeteorEndNarrow-possibleMeteorStartNarrow;
            %see if the length of this meteor is better than the best
            if possibleMeteorLengthNarrow > possibleMeteorLengthBestNarrow
                %update the best possible meteor
                possibleMeteorStartBestNarrow = possibleMeteorStartNarrow;
                possibleMeteorEndBestNarrow = possibleMeteorEndNarrow;
                possibleMeteorLengthBestNarrow = possibleMeteorLengthNarrow;
            end
        end
    end
end

if possibleMeteorLengthBestNarrow > possibleMeteorLengthBest
    disp('Narrow Beam has more points than the Wide Beam. Weird!')
    possibleMeteorStartBest = possibleMeteorStartBestNarrow;
    possibleMeteorEndBest = possibleMeteorEndBestNarrow;
    possibleMeteorLengthBest = possibleMeteorLengthBestNarrow;
    radarPowerForRange = narrowBeamPower;
    timeForRange = narrowBeamTime;
else
    %the wide beam is the correct radar beam to use for range extrapolation
    radarPowerForRange = wideBeamPower;
    timeForRange = wideBeamTime;
end
% Now we have the indeces for the meteor in the radar data, find the range
% information
radarGoodIndeces = linspace(possibleMeteorStartBest,possibleMeteorEndBest,possibleMeteorLengthBest+1);
radarGoodRanges = zeros(1,numel(radarGoodIndeces));
radarGoodPowers = zeros(1,numel(radarGoodIndeces));
for i = 1:numel(radarGoodIndeces)
   %get the index of maximum power at the good indeces
   [~,maxPowerIndex] = max(radarPowerForRange(:,radarGoodIndeces(i)));
   %get the range value at the index
   radarGoodRanges(i) = wideBeamRange(maxPowerIndex);
   radarGoodPowers(i) = radarPowerForRange(maxPowerIndex,radarGoodIndeces(i));
end
%find the indeces of the times we care about
radarValidTimeIndeces = zeros(1,numel(meteorTimes));
meteorRange = zeros(1,numel(meteorTimes));
for i = 1:numel(meteorTimes)
   [~,radarValidTimeIndeces(i)] = min(abs(timeForRange-meteorTimes(i)));
end
%get an interpolation for the data for the frames we care about
%meteorTimeFit = polyfit(wideBeamTime(radarGoodIndeces),radarGoodRanges,1);
%meteorRanges = meteorTimeFit(1)*meteorTimes+meteorTimeFit(2);
meteorTimeFit = polyfit(timeForRange(radarGoodIndeces),radarGoodRanges,pRangeFit);
meteorRanges = polyval(meteorTimeFit,meteorTimes);

%(TESTING) Plot the meteor ranges from the radar, and the extrapolated
%values on the same graph
figure(4)
movegui('southwest')
minTime = min([meteorTimes(1) wideBeamTime(radarGoodIndeces(1))]);
plot((meteorTimes-minTime)*matlabTmsInSec,meteorRanges,'r')
hold on
plot((wideBeamTime(radarGoodIndeces)-minTime)*matlabTmsInSec,radarGoodRanges)
legend('Optical','Radar')
xlabel('Time (sec)')
ylabel('Range (km)')
%hold off
title('Meteor Ranges From Radar and Extrapolated to Optical Time')
%(END TESTING)

%so now we have the range, azimuth, and elevation of the meteor
%let's find the velocity!
%convert to cartesian coordinates
meteorX = meteorRanges.*sind(90-meteorEl).*cosd(meteorAz);
meteorY = meteorRanges.*sind(90-meteorEl).*sind(meteorAz);
meteorZ = meteorRanges.*cosd(90-meteorEl);
meteorVel = sqrt((meteorX(2:end)-meteorX(1:end-1)).^2+(meteorY(2:end)-meteorY(1:end-1)).^2+(meteorZ(2:end)-meteorZ(1:end-1)).^2)./(meteorTimes(2:end)-meteorTimes(1:end-1))./matlabTmsInSec;
%(TESTING) Plot the meteor position
figure(5)
movegui('south')
plot3(meteorX,meteorY,meteorZ)
hold on
plot3(meteorX,meteorY,meteorZ,'x')
for i = 1:numel(meteorX)
    text(meteorX(i),meteorY(i),meteorZ(i), num2str(meteorFrame(i)), 'VerticalAlignment','bottom', ...
                                                            'HorizontalAlignment','right')
    text(meteorX(i),meteorY(i),meteorZ(i), ['\color{red} ' num2str(round(meteorWideGain(i)))],'VerticalAlignment','top', ...
                                                            'HorizontalAlignment','left')                                      
end
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Altitude (km)')
grid on
%hold off
title('Meteor Position')
%plot the meteor velocity vs time
figure(6)
movegui('southeast')
plot((meteorTimes(1:end-1)-meteorTimes(1))*matlabTmsInSec,meteorVel)
ylabel('Velocity (km/sec)')
xlabel('Time (sec)')
title('Meteor Velocity')
%(END TESTING)
mean(meteorVel)
std(meteorVel)

%plot the radar data
figure(7)
movegui('northwest')
subplot(2,1,1)
imagesc((wideBeamTime-wideBeamTime(1))*matlabTmsInSec,wideBeamRange,db(wideBeamPower))
xlabel('Time (sec)')
ylabel('Range (km)')
title('Wide Beam')
set(gca,'ydir','normal')
subplot(2,1,2)
imagesc((narrowBeamTime-wideBeamTime(1))*matlabTmsInSec,narrowBeamRange,db(narrowBeamPower))
set(gca,'ydir','normal')
xlabel('Time (sec)')
ylabel('Range (km)')
title('Narrow Beam')


% END OF NEW METHOD 3/9/2015

% Plot movie
if show_SBR == 1
    figure()
    for i = 1:size(SBR,3)
        imagesc(SBR(:,:,i))
        set(gca,'clim',[0 max(SBR(:))])
        hold on
        plot(xCoords(i),yCoords(i),'rx')
        hold off
        title('Signal to Background Ratio (dB)')
        xlabel(datestr(tUTCU(i),'dd-mmm-yyyy HH:MM:SS.FFF'))
        colormap(gray)
        colorbar
        pause(plotPauseTime)
    end
end

if show_signal == 1
    figure()
    for i = 1:size(SBR,3)
        imagesc(signalUltra(:,:,i))
        hold on
        plot(xCoords(i),yCoords(i),'rx')
        hold off
        title('Signal')
        xlabel(datestr(tUTCU(i),'dd-mmm-yyyy HH:MM:SS.FFF'))
        colormap(gray)
        set(gca,'clim',[0 max(signalUltra(:))/10])
        colorbar
        pause(plotPauseTime)
    end
end

% for i = 1:size(SBR,3)
%         imagesc(dataUltra(:,:,i))
%         title('Signal')
%         xlabel(datestr(tUTCU(i),'dd-mmm-yyyy HH:MM:SS.FFF'))
%         colorbar
%         pause(0.1)
% end

%% Get the RCS of the head echo
powerReceived = radarGoodPowers;
rcs = powerReceived*4*pi*range.^2/(radarGain*powerTransmit);
end