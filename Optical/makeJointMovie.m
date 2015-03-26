
% %make a movie
% %OPTICAL MOVIE
% fig1 = figure(1);
% winsize = get(fig1,'Position');
% winsize(1:2) = [0 0];
% numframesUltra = stopFrame-startFrame+1;
% movieUltra = moviein(numframesUltra,fig1,winsize);
% set(fig1,'NextPlot','replacechildren')
% for i = startFrame:stopFrame
%     imagesc(10*log10(double(dataUltra(:,:,i))./dataUltraBackground),[minSNROptical maxSNROptical])
%     colormap(gray)
%     title(['Frame' num2str(i)])
%     colorbar()
%     movieUltra(:,i-startFrame+1) = getframe(fig1,winsize);
% end
% 
% %RADAR MOVIE
% figRadar = figure();
% winsize = get(figRadar,'Position');
% winsize(1:2) = [0 0];
% numframes = size(opticalToRadarFrameIndeces,2)+1;
% movieRadar = moviein(numframes,figRadar,winsize);
% set(figRadar,'NextPlot','replacechildren')
% %first plot a blank image with the correct size
% subplot(2,1,1)
% radarBuiltUp = ones(2,size(rng2,1),stopRadarFrameIndex-startRadarFrameIndex+1);
% radarBuiltUp(1,:,:) = radarBuiltUp(1,:,:)*minSNRRadar(1);
% radarBuiltUp(2,:,:) = radarBuiltUp(2,:,:)*minSNRRadar(2);
% imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),...
%     rng2/1e3,squeeze(radarBuiltUp(1,:,:)),[minSNRRadar(1),maxSNRRadar(1)])
% colorbar()
% subplot(2,1,2)
% imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),...
%     rng2/1e3, squeeze(radarBuiltUp(2,:,:)),[minSNRRadar(2),maxSNRRadar(2)])
% colorbar()
% noiseRadar = [median(median(squeeze(odat(1,:,:)))), median(median(squeeze(odat(2,:,:))))];
% for i = 1:size(opticalToRadarFrameIndeces,2)
%     startIndex = opticalToRadarFrameIndeces(1,i);
%     stopIndex = opticalToRadarFrameIndeces(2,i);
%     
%     %plot the strip of radar data
%     for ibm = 1:2
%         radarBuiltUp(ibm,:,startIndex-startRadarFrameIndex+1:stopIndex-startRadarFrameIndex+1) = ...
%             10*log10(squeeze(odat(ibm,:,startIndex:stopIndex))./noiseRadar(ibm));
%         subplot(2,1,ibm);
%         imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),rng2/1e3,...
%             squeeze(radarBuiltUp(ibm,:,:)),[minSNRRadar(ibm),maxSNRRadar(ibm)]);
%         colorbar()
%         set(gca,'ydir','normal');
%         %ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
%         ylabel(['Range (km);' ibms(ibm) ' Beam'])
%     end
%     movieRadar(:,i) = getframe(figRadar,winsize);
% end


figJoint = figure('units','pixels','outerposition',[0 0 900 900]);
winsize = get(figJoint,'Position');
winsize(1:2) = [0 0];
numframesJoint = stopFrame-startFrame+1;
movieJoint = moviein(numframesJoint,figJoint,winsize);
set(figJoint,'NextPlot','replacechildren')
radarBuiltUp = ones(2,size(rng2,1),stopRadarFrameIndex-startRadarFrameIndex+1);
radarBuiltUp(1,:,:) = radarBuiltUp(1,:,:)*minSNRRadar(1);
radarBuiltUp(2,:,:) = radarBuiltUp(2,:,:)*minSNRRadar(2);
subplot(4,2,[5 6])
imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),...
    rng2/1e3,squeeze(radarBuiltUp(1,:,:)),[minSNRRadar(1),maxSNRRadar(1)])
colorbar()
subplot(4,2,[7 8])
imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),...
    rng2/1e3, squeeze(radarBuiltUp(2,:,:)),[minSNRRadar(2),maxSNRRadar(2)])
colorbar()
noiseRadar = [median(median(squeeze(odat(1,:,:)))), median(median(squeeze(odat(2,:,:))))];
for i = startFrame:stopFrame
    frameIndex = i;
    startIndex = opticalToRadarFrameIndeces(1,i-startFrame+1);
    stopIndex = opticalToRadarFrameIndeces(2,i-startFrame+1);
    %plot the ultra data
    subplot(4,2,[1 2 3 4])
    imagesc(10*log10(double(dataUltra(:,:,frameIndex))./dataUltraBackground),[minSNROptical maxSNROptical])
    %imagesc(10*log10(double(dataUltra(:,:,frameIndex))./dataUltraBackground))
    colormap(gray)
    freezeColors
    title(strcat(datestr(tUTCU(i),'mmmm dd, yyyy HH:MM:SS.FFF'), ' UTC'))
    %title(strcat('Frame ', num2str(i)))
    colorbar()
    cbfreeze()
    %plot the radar data
    for ibm = 1:2
        radarBuiltUp(ibm,:,startIndex-startRadarFrameIndex+1:stopIndex-startRadarFrameIndex+1) = ...
            10*log10(squeeze(odat(ibm,:,startIndex:stopIndex))./noiseRadar(ibm));
        subplot(4,2,[ibm*2+3 ibm*2+4]);
        imagesc(tms(1,startRadarFrameIndex:stopRadarFrameIndex)-tms(1,startRadarFrameIndex),rng2/1e3,...
            squeeze(radarBuiltUp(ibm,:,:)),[minSNRRadar(ibm),maxSNRRadar(ibm)]);
        colormap(jet)
        freezeColors
        colorbar()
        cbfreeze()
        set(gca,'ydir','normal');
        ylabel('Range (km)')
        title(strcat(ibms(ibm), ' Beam'))
        if ibm == 2
           xlabel('Time (sec)') 
        end
    end
    
    if i == 37
        pause(0.2)
    end
    if i == 42
        pause(0.2)
    end
    movieJoint(:,i-startFrame+1) = getframe(figJoint,winsize);
end

%movie2avi(movieUltra, 'ultraMovie.avi', 'FPS', 53)
%movie2avi(movieRadar, 'radarMovie.avi', 'FPS', 53)
movie2avi(movieJoint, strcat(filenameOptical(1:end-9),'compressed','.avi'), 'FPS', 26, 'compression', 'FFDS')
%movie2avi(movieJoint, strcat(filenameOptical(1:end-9),'uncompressed','.avi'), 'FPS', 26)
%movie(figRadar,movieRadar,1,53,winsize)
%movie(figJoint,movieJoint,1,53,winsize)
movie(figJoint,movieJoint,1,26,winsize)