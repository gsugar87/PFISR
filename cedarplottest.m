close all

load 'testBoth.mat'

titleFont = 20;
axisFont = 16;
tickFont = 15;

%imagesc(fliplr(10*log10(sum_data_streak_ixon'./median(median(sum_data_streak_ixon)))))
xu = flipud(sum_data_streak_ultra);
imagesc(10*log10(xu./median(median(sum_data_streak_ultra))))
set(gca,'ydir','normal');
set(gca,'FontSize',tickFont)
title('No Filter','fontsize',titleFont)
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','SNR (dB)','fontsize',axisFont)
set(cbar_handle,'FontSize',tickFont)
caxis([0,3])
%headCoordinates = ginput(2);
headCoordinates = [[171, 392];[243,370]]
slope = (headCoordinates(2,2)-headCoordinates(1,2))/(headCoordinates(2,1)-headCoordinates(1,1));
xStart = round(headCoordinates(1,1));
xEnd = round(headCoordinates(2,1));
powerUltra = zeros(xEnd-xStart+1,1);

for i = xStart:xEnd
    yIndex = round(headCoordinates(1,2)+(i-xStart)*slope);
    powerUltra(i-xStart+1) = max(xu(yIndex-2:yIndex+2,i));
end

xi = flipud(fliplr(sum_data_streak_ixon'));
figure;
imagesc(10*log10(xi./median(median(sum_data_streak_ixon))))
set(gca,'ydir','normal');
set(gca,'FontSize',tickFont)
title('> 475nm Filter','fontsize',titleFont)
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','SNR (dB)','fontsize',axisFont)
set(cbar_handle,'FontSize',tickFont)
caxis([0,0.03])
%headCoordinates = ginput(2);
%headCoordinates = [[205,395],[240,381]]
headCoordinates = [[186,402];[258,374]];
slope = (headCoordinates(2,2)-headCoordinates(1,2))/(headCoordinates(2,1)-headCoordinates(1,1));
xStart = round(headCoordinates(1,1));
xEnd = round(headCoordinates(2,1));
powerIxon = zeros(xEnd-xStart+1,1);

for i = xStart:xEnd
    yIndex = round(headCoordinates(1,2)+(i-xStart)*slope);
    powerIxon(i-xStart+1) = max(xi(i,yIndex-2:yIndex+2));
end

%find the maximum index of each power
[dummy, maxIxonIndex] = max(powerIxon);
[dummy, maxUltraIndex] = max(powerUltra);
figure
titleFont = 20;
axisFont = 16;
tickFont = 15;
plot(powerUltra./powerIxon)
set(gca,'FontSize',tickFont)
xlabel('Pixel','fontsize',axisFont);
ylabel('No Filter/Filter Intensity','fontsize',axisFont);
title('No Filter/Filter Camera Ratio')


