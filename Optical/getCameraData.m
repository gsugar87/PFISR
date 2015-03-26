%% getCameraData
% This script will read in an image from either the Ixon Classic camera or
% the Ixon Ultra camera and plot the data
close all
%fileName = 'E:\PFISR Images\UltraPFRR\2014-03-30\2014-03-30T10-46-CamSer7196.DMCdata';
%fileName = 'G:\PFISR Images\UltraPFRR\2014-03-30\2014-03-30T10-46-CamSer7196.DMCdata';
fileNameIxon = 'E:\PFISR Images\Ixon\2014-03-30\2014-03-30T10-58-CamSer1387.DMCdata';
fileNameUltra = 'E:\PFISR Images\UltraPFRR\2014-03-30\2014-03-30T10-46-CamSer7196.DMCdata';
FPSIxon = 33.00125;
FPSUltra = 53.00125;
secondsToCheck = 2;
secondsOfBackground = 1.5;
%desiredTime = [3,30,13,11,43];
%desiredTime = [3,30,12,39,03];
%desiredTime = [3,30,11,24,38.5]; %this one in both!
%desiredTime = [3,30,11,48,21];
%desiredTime = [3,30,11,54,00];
desiredTime = [3,30,10,49,38];
desiredTime = desiredTime(1)*30*24*60*60+desiredTime(2)*24*60*60+desiredTime(3)*60*60+desiredTime(4)*60+desiredTime(5);
%frameStart = 1;
%frameEnd = 3;
frameStartI = 1000;
frameStartU = 1000;
frameEndI = ceil(frameStartI + FPSIxon);
frameEndU = ceil(frameStartU + FPSUltra);
xPix = 512;
yPix = 512;
xBin = 1;
yBin = 1;
%playMovie = 0.01;
playMovie = 0;
ClimU = [100,1100];
ClimI = [960,1070];
%read in the calibration files
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\ixonCalibrate2.mat')
azI = rot90(az,3);
elI = rot90(el,3);
load('C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\ultraCalibrate2.mat')
azU = az;
elU = el;
% rawDMCreader(BigFN,xPix,yPix,xBin,yBin,FrameInd,playMovie,Clim,rawFrameRate,startUTC)
[dataI,~,tUTCI] = rawDMCreaderGlenn(fileNameIxon,xPix,yPix,xBin,yBin,frameStartI:frameEndI,0,ClimI,'auto','auto');
[dataU,~,tUTCU] = rawDMCreaderGlenn(fileNameUltra,xPix,yPix,xBin,yBin,frameStartU:frameEndU,0,ClimU,'auto','auto');
% Next lines are for saving the calibrated image, which was already
% created on 7/23/2014 by Glenn
% dataSum = sum(dataIxonBackground,3);
% dataImage = log10(dataSum);
% dataImage = dataImage - min(min(dataImage));
% dataImage = dataImage./max(max(dataImage));
% imwrite(dataImage, 'IxonCalibrate.png', 'png')

frameOneTimeI = datestr(tUTCI(1));
frameOneSecI = str2num(frameOneTimeI(end-1:end));
frameOneMinI = str2num(frameOneTimeI(end-4:end-3));
frameOneHrI = str2num(frameOneTimeI(end-7:end-6));
frameOneDayI = str2num(frameOneTimeI(1:2));
frameOneMonthI = frameOneTimeI(4:6);
if frameOneMonthI == 'Mar'
    frameOneMonthI = 3.0;
else
    frameOneMonthI = 4.0;
end
frameOneTimeI = frameOneMonthI*30*24*60*60 + frameOneDayI*24*60*60 + frameOneHrI*60*60 + ...
    frameOneMinI*60 + frameOneSecI;
frameOneTimeU = datestr(tUTCU(1));
frameOneSecU = str2num(frameOneTimeU(end-1:end));
frameOneMinU = str2num(frameOneTimeU(end-4:end-3));
frameOneHrU = str2num(frameOneTimeU(end-7:end-6));
frameOneDayU = str2num(frameOneTimeU(1:2));
frameOneMonthU = frameOneTimeU(4:6);
if frameOneMonthU == 'Mar'
    frameOneMonthU = 3.0;
else
    frameOneMonthU = 4.0;
end
frameOneTimeU = frameOneMonthU*30*24*60*60 + frameOneDayU*24*60*60 + frameOneHrU*60*60 + ...
    frameOneMinU*60 + frameOneSecU;

%get number of seconds between the first frame and the desired frame
secondsIntoFileI = desiredTime-frameOneTimeI+frameStartI/FPSIxon;
secondsIntoFileU = desiredTime-frameOneTimeU+frameStartU/FPSUltra;
if secondsIntoFileI > 0
    framesIntoFileI = ceil(secondsIntoFileI*FPSIxon);
else
    framesIntoFileI = 0;
    'Check the desired time'
end
if secondsIntoFileU > 0
    framesIntoFileU = ceil(secondsIntoFileU*FPSUltra);
else
    framesIntoFileU = 0;
    'Check the desired time'
end

%% MANUAL SET FRAMES INTO FILE!!! %%%%
disp('THE CHOSEN FRAMES ARE AT LINES 97-103!!!')
frameBuffer = 20;
%fall = [31579 35833];
%saveFilename = '0330105633ultra.mat';
fall = [357521 357546];
saveFilename = '0330123903ultra.mat';
finit = fall(1);
fend = fall(2);
framesIntoFileU = finit-frameBuffer;
[dut,~,tUTCUtt] = rawDMCreaderGlenn(fileNameUltra,xPix,yPix,xBin,yBin,finit:fend,1,ClimU,'auto','auto');

%'Getting the Classic Data'
%[dataIxon,~,tUTCI] = rawDMCreaderGlenn(fileNameIxon,xPix,yPix,xBin,yBin,framesIntoFileI:framesIntoFileI+ceil(secondsToCheck*FPSIxon),playMovie,ClimI,'auto','auto');
%'Getting the Classic Background'
%[dataIxonBackground,~,tUTCIBackground] = rawDMCreaderGlenn(fileNameIxon,xPix,yPix,xBin,yBin,framesIntoFileI-ceil(FPSIxon*secondsOfBackground)-1:framesIntoFileI-1,playMovie,ClimI,'auto','auto');
%backgroundStdIxonRot = rot90(std(double(dataIxonBackground),[],3));
%dataIxonBackgroundRot = rot90(mean(dataIxonBackground,3));
'Getting the Ultra Data'
%[dataUltra,~,tUTCU] = rawDMCreaderGlenn(fileNameUltra,xPix,yPix,xBin,yBin,framesIntoFileU:framesIntoFileEndU,playMovie,ClimU,'auto','auto');
[dataUltra,~,tUTCU] = rawDMCreaderGlenn(fileNameUltra,xPix,yPix,xBin,yBin,framesIntoFileU:framesIntoFileU+secondsToCheck*FPSUltra,playMovie,ClimU,'auto','auto');
'Getting the Ultra Background'
[dataUltraBackground,~,tUTCUBackground] = rawDMCreaderGlenn(fileNameUltra,xPix,yPix,xBin,yBin,framesIntoFileU-ceil(FPSUltra*secondsOfBackground)-1:framesIntoFileU-1,playMovie,ClimU,'auto','auto');
backgroundStdUltra = std(double(dataUltraBackground),[],3);
dataUltraBackground = mean(dataUltraBackground,3);
dataMovie = dataUltra;
dataMovieBackground = dataUltraBackground;
framesIntoMovieFile = framesIntoFileU;
% dataMovie = dataIxon;
% dataMovieBackground = dataIxonBackground;
% framesIntoMovieFile = framesIntoFileI;
% for i = 1:size(dataMovie,3)
%     imagesc(double(dataMovie(:,:,i)) - dataMovieBackground)
%     frameNum = (i-1)+framesIntoMovieFile;
%     title(['Ultra or Ixon Frame: ' num2str(frameNum)])
%     colorbar
%     pause(0.3)
% end

%get various images
%sum
dataSumU = sum(dataUltra,3);
dataStdU = std(double(dataUltra),[],3);
meteorImageU = dataStdU./mean(dataUltra,3);
%dataSumIRot = rot90(sum(dataIxon,3),3);
%dataStdI = std(double(dataIxon),[],3);
%meteorImageIRot = rot90(dataStdI./mean(dataIxon,3),3);

%get the coordinates of the meteor box
figure()
imagesc(meteorImageU)
title('Ultra Standard Deviation / Mean Intensity')
colormap(gray)
% figure()
% colormap(gray)
% h = imagesc(log10(rot90(dataSumI,3)));
%get coordinates of the meteor box
[xU,yU,button] = ginput(2);
if button(1) ~= 1 || button(2) ~= 1
    %we want to zoom in on the image
    xU = round(xU);
    yU = round(yU);
    xU(xU < 1) = 1;
    xU(xU > size(meteorImageU,1)) = size(meteorImageU,1);
    yU(yU < 1) = 1;
    yU(yU > size(meteorImageU,1)) = size(meteorImageU,1);
    xU = [min(xU),max(xU),max(xU),min(xU)];
    yU = [min(yU),min(yU),max(yU),max(yU)];
    meteorImageUZoom = meteorImageU(min(yU):max(yU),min(xU):max(xU));
    figure()
    imagesc(meteorImageUZoom)
    title('Ultra Standard Deviation / Mean Intensity')
    colormap(gray)
    [xUZoom,yUZoom] = ginput(2);
    xUZoom = round(xUZoom);
    yUZoom = round(yUZoom);
    xUZoom(xUZoom < 1) = 1;
    xUZoom(xUZoom > size(meteorImageUZoom,2)) = size(meteorImageUZoom,2);
    yUZoom(yUZoom < 1) = 1;
    yUZoom(yUZoom > size(meteorImageUZoom,1)) = size(meteorImageUZoom,1);
    xU = [xUZoom(1)+xU(1)-1,xUZoom(2)+xU(1)-1,xUZoom(2)+xU(1)-1,xUZoom(1)+xU(1)-1];
    yU = [yUZoom(1)+yU(1)-1,yUZoom(1)+yU(1)-1,yUZoom(2)+yU(1)-1,yUZoom(2)+yU(1)-1];
else
    
    xU = round(xU);
    yU = round(yU);
    xU(xU < 1) = 1;
    xU(xU > size(meteorImageU,1)) = size(meteorImageU,1);
    yU(yU < 1) = 1;
    yU(yU > size(meteorImageU,1)) = size(meteorImageU,1);
    xU = [min(xU),max(xU),max(xU),min(xU)];
    yU = [min(yU),min(yU),max(yU),max(yU)];
end
%get the total intensity over both Ixon and Ultra meteor rectangles
%go through each column and find where the intensity falls below the noise
meteorIntensityUltra2D = dataSumU(min(yU):max(yU),min(xU):max(xU))- ...
    dataUltraBackground(min(yU):max(yU),min(xU):max(xU)).*size(dataUltra,3);
%%%%%%%%%%% ULTRA INTENSITY START %%%%%%%%%%%%%%%%%%
meteorIntensityUltra = 0;
meteorSignalUltraMax = 0;
meteorIntensityUltraError = 0;
slopeUltra = (max(xU)-min(xU))/(max(yU)-min(yU));
%see if the meteor is going diagonal up or down
%look at top left and compare signal strength to bottom left
diagDown = zeros(max(yU)-min(yU),1);
diagUp = zeros(max(yU)-min(yU),1);
for rowIndex = 1:size(diagDown,1)
    diagDown(rowIndex) = meteorImageU(min(yU)+rowIndex-1,min(xU) + round((rowIndex-1)*slopeUltra));
    diagUp(rowIndex) = meteorImageU(min(yU)+rowIndex-1,max(xU) - round((rowIndex-1)*slopeUltra));
end
if median(diagDown) > median(diagUp)
    meteorDirection = 1; %the meteor is diagonal down, the usual
else
    meteorDirection = -1;
end
%go through each row of the image
for rowIndex = 1:size(meteorIntensityUltra2D,1)
    %look through each row index and find the center of the meteor for that
    %row
    if meteorDirection == 1
        xIndex = 1+round((rowIndex-1)*slopeUltra);
    else
        xIndex = size(meteorIntensityUltra2D,2) - round((rowIndex-1)*slopeUltra);
    end
    intensityRow = meteorIntensityUltra2D(rowIndex,:);
    %look above xIndex
    xIndexUp = xIndex;
    while (intensityRow(xIndexUp) > 0) && (xIndexUp < size(intensityRow,2))
        %increase the up index until it can't be increased any more
        xIndexUp = xIndexUp + 1;
    end
    %look below xIndex
    xIndexDown = xIndex;
    while intensityRow(xIndexDown) > 0 && xIndexDown > 1
        xIndexDown = xIndexDown - 1;
    end
    meteorIntensityUltra = meteorIntensityUltra + sum(intensityRow(xIndexDown:xIndexUp));
    meteorSignalUltraMax = max([intensityRow(xIndexDown:xIndexUp),meteorSignalUltraMax]);
    meteorIntensityUltraError = sqrt(meteorIntensityUltraError^2 + sum(intensityRow(xIndexDown:xIndexUp).^2));
end
disp('Ultra max:')
disp(meteorSignalUltraMax)
disp(meteorIntensityUltra)
disp(meteorIntensityUltraError)
%%%%%%%%%%% ULTRA INTENSITY FINISHED %%%%%%%%%%%%%%%%%%

%get az and el limits
az = [azU(xU(1),yU(1)),azU(xU(2),yU(2)),azU(xU(3),yU(3)),azU(xU(4),yU(4))];
el = [elU(xU(1),yU(1)),elU(xU(2),yU(2)),elU(xU(3),yU(3)),elU(xU(4),yU(4))];
xI = zeros(4,1);
yI = zeros(4,1);

% % do the same for the ixon camera
% figure()
% imagesc(meteorImageIRot)
% title('Classic (rotated) Standard Deviation / Mean Intensity')
% colormap(gray)
% [xI,yI,button] = ginput(2);
% if button(1) == 3 || button(2) == 3
%     %we want to zoom in on the image
%     xI = round(xI);
%     yI = round(yI);
%     xI(xI < 1) = 1;
%     xI(xI > size(meteorImageIRot,1)) = size(meteorImageIRot,1);
%     yI(yI < 1) = 1;
%     yI(yI > size(meteorImageIRot,1)) = size(meteorImageIRot,1);
%     xI = [min(xI),max(xI),max(xI),min(xI)];
%     yI = [min(yI),min(yI),max(yI),max(yI)];
%     meteorImageIRotZoom = meteorImageIRot(min(yI):max(yI),min(xI):max(xI));
%     figure()
%     imagesc(meteorImageIRotZoom)
%     title('Classic (rotated) Standard Deviation / Mean Intensity')
%     colormap(gray)
%     [xIZoom,yIZoom] = ginput(2);
%     xIZoom = round(xIZoom);
%     yIZoom = round(yIZoom);
%     xIZoom(xIZoom < 1) = 1;
%     xIZoom(xIZoom > size(meteorImageIRotZoom,2)) = size(meteorImageIRotZoom,2);
%     yIZoom(yIZoom < 1) = 1;
%     yIZoom(yIZoom > size(meteorImageIRotZoom,1)) = size(meteorImageIRotZoom,1);
%     xI = [xIZoom(1)+xI(1)-1,xIZoom(2)+xI(1)-1,xIZoom(2)+xI(1)-1,xIZoom(1)+xI(1)-1];
%     yI = [yIZoom(1)+yI(1)-1,yIZoom(1)+yI(1)-1,yIZoom(2)+yI(1)-1,yIZoom(2)+yI(1)-1];
% elseif button(1) == 1 || button(2) == 1
%     %get az and el limits
%     xI = round(xI);
%     yI = round(yI);
%     xI(xI < 1) = 1;
%     xI(xI > size(meteorImageIRot,1)) = size(meteorImageIRot,1);
%     yI(yI < 1) = 1;
%     yI(yI > size(meteorImageIRot,1)) = size(meteorImageIRot,1);
%     xI = [min(xI),max(xI),max(xI),min(xI)];
%     yI = [min(yI),min(yI),max(yI),max(yI)];
% else
%     %this is activated when a middle mouse button was pressed
%     %we can't see the meteor in the ixon camera, so use the az/el values
%     %from the ultra data
%     %find the corresponding points for the ixon camera using 
%     for point = 1:4
%         minError = 100; %initialize a large min error
%         for xIndex = 1:size(azI,1)
%             for yIndex = 1:size(elI,1)
%                 if ((azI(xIndex,yIndex)-az(point))^2 + (elI(xIndex,yIndex)-el(point))^2) < minError
%                     minError = (azI(xIndex,yIndex)-az(point))^2 + (elI(xIndex,yIndex)-el(point))^2;
%                     xI(point) = xIndex;
%                     yI(point) = yIndex;
%                 end
%             end
%         end
%     end
% end
%        
% 
% %get the total intensity over both Ixon and Ultra meteor rectangles
% %go through each column and find where the intensity falls below the noise
% meteorIntensityIxon2D = dataSumIRot(min(yI):max(yI),min(xI):max(xI))- dataIxonBackgroundRot(min(yI):max(yI),min(xI):max(xI)).*size(dataIxon,3);
% %%%%%%%%%%% IXON INTENSITY START %%%%%%%%%%%%%%%%%%
% meteorIntensityIxon = 0;
% meteorSignalIxonMax = 0;
% meteorIntensityIxonError = 0;
% slopeIxon = (max(xI)-min(xI))/(max(yI)-min(yI));
% for rowIndex = 1:size(meteorIntensityIxon2D,1)
%     %look through each row index
%     xIndex = 1+round((rowIndex-1)*slopeIxon);
%     intensityRow = meteorIntensityIxon2D(rowIndex,:);
%     %look above xIndex
%     xIndexUp = xIndex;
%     while (intensityRow(xIndexUp) > 0) && (xIndexUp < size(intensityRow,2))
%         %increase the up index until it can't be increased any more
%         xIndexUp = xIndexUp + 1;
%     end
%     %look below xIndex
%     xIndexDown = xIndex;
%     while intensityRow(xIndexDown) > 0 && xIndexDown > 1
%         xIndexDown = xIndexDown - 1;
%     end
%     meteorIntensityIxon = meteorIntensityIxon + sum(intensityRow(xIndexDown:xIndexUp));
%     meteorSignalIxonMax = max([intensityRow(xIndexDown:xIndexUp),meteorSignalIxonMax]);
%     meteorIntensityIxonError = sqrt(meteorIntensityIxonError^2 + sum(intensityRow(xIndexDown:xIndexUp).^2));
% end
% disp('Ixon max')
% disp(meteorSignalIxonMax)
% disp(meteorIntensityIxon)
% disp(meteorIntensityIxonError)
% %%%%%%%%%%% IXON INTENSITY FINISHED %%%%%%%%%%%%%%%%%%
% meteorIntensityRatio = meteorIntensityIxon/meteorIntensityUltra;
% meteorIntensityRatioError = abs(meteorIntensityRatio)*sqrt((meteorIntensityIxon/meteorIntensityIxonError)^2 + (meteorIntensityUltra/meteorIntensityUltraError)^2);
% disp(sprintf('The meteor intensity ratio is: %f +/- %f', [meteorIntensityRatio, meteorIntensityRatioError]));

intensityUltra = sum(sum(dataSumU(min(yU):max(yU),min(xU):max(xU))- dataUltraBackground(min(yU):max(yU),min(xU):max(xU)).*size(dataUltra,3)));
intensityUltraError = sqrt(sum(sum(backgroundStdUltra.^2))*size(dataUltra,3));
% intensityIxon = sum(sum(dataSumIRot(min(yI):max(yI),min(xI):max(xI))- dataIxonBackgroundRot(min(yI):max(yI),min(xI):max(xI)).*size(dataIxon,3)));
% intensityIxonError = sqrt(sum(sum(backgroundStdIxonRot(min(yI):max(yI),min(xI):max(xI)).^2))*size(dataIxon,3));
% intensityRatio = intensityIxon/intensityUltra;
% intensityRatioError = abs(intensityRatio)*sqrt((intensityIxonError/intensityIxon)^2+(intensityUltraError/intensityUltra)^2);
% disp(sprintf('The intensity ratio is: %f +/- %f', [intensityRatio, intensityRatioError]));
% figure()
% imagesc(meteorImageIRot)
% title('Classic (rotated) Standard Deviation / Mean Intensity')
% colormap(gray)
% figure()
% imagesc(meteorImageU(min(yU):max(yU),min(xU):max(xU)))
% title('Ultra Standard Deviation / Mean Intensity Zoom')
% colormap(gray)
% figure()
% imagesc(meteorImageIRot(min(yI):max(yI),min(xI):max(xI)))
% title('Classic (rotated) Standard Deviation / Mean Intensity Zoom')
% colormap(gray)
% figure()
% imagesc(sum(dataSumIRot(min(yI):max(yI),min(xI):max(xI),:),3))
% colorbar()
% title('Classic (rotated) Total Intensity Zoom')
% colormap(gray)
figure()
imagesc(sum(dataSumU(min(yU):max(yU),min(xU):max(xU),:),3))
colorbar()
title('Ultra Total Intensity Zoom')
colormap(gray)
% figure()
% imagesc(meteorIntensityIxon2D,[0 max(max(meteorIntensityIxon2D))])
% colorbar()
% title('Ixon (rotated) Total Signal Zoom')
% colormap(gray)
% colorbar()
figure()
imagesc(meteorIntensityUltra2D,[0 max(max(meteorIntensityUltra2D))])
colorbar()
title('Ultra Total Signal Zoom')
colormap(gray)
colorbar()
figure()
imagesc(10*log10(dataSumU(min(yU):max(yU),min(xU):max(xU))./ ...
    (dataUltraBackground(min(yU):max(yU),min(xU):max(xU)).*size(dataUltra,3))))
title('Ultra SNR Zoom')
colorbar()
% figure()
% imagesc(10*log10(dataSumIRot(min(yI):max(yI),min(xI):max(xI))./ ...
%     (dataIxonBackgroundRot(min(yI):max(yI),min(xI):max(xI)).*size(dataIxon,3))))
% title('Classic (rotated) SNR Zoom')
% colorbar()

%save image as png
if 1==0
    dataImage = log10(dataSum);
    dataImage = dataImage - min(min(dataImage));
    dataImage = dataImage./max(max(dataImage));
    imwrite(dataImage, 'test.png', 'png')
end

%save the data
if 1 == 1
   save(['C:\Users\Glenn\Documents\MATLAB\PFISR\Optical\Meteors\' saveFilename], 'dataUltra','dataUltraBackground','tUTCU')
end