close all
time = 1396170000;
time = time - 9*60*60; %march 30 2014
month = 3;
day = 31;
%time = time + 24*60*60; %march 31 2014 UTC
% hour = 8;
% minute = 25;
% second = 52;
% hour = 11;
% minute = 20;
% second = 52;
% hour = 11;
% minute = 29;
% second = 27;
%hour = 11;
%minute = 00;
%second = 16;
%eventTime = [13, 32, 43];
%eventTime = [11, 00, 15];
%eventTime = [11, 16, 37];
%eventTime = [11, 24, 40];
%eventTime = [11, 48, 22];
%eventTime = [12, 13, 07];
%eventTime = [12, 23, 29];
%eventTime = [12, 39, 04];
%eventTime = [10, 56, 34];
eventTime = [10,57,54];
hour = eventTime(1);
minute = eventTime(2);
second = eventTime(3);
time = time + hour*60*60+minute*60+second;

meteor_name = ['0' num2str(month) num2str(day) num2str(hour) num2str(minute) num2str(second)];

timezone = -8;
%ddir = '/Volumes/AMISR_019/Data AMISR Poker/20140327.017/';
%ddir = 'F:\20140331.005\';
%ddir = 'C:\PFISR\Data\';
%ddir = 'E:\20140330.007\';
ddir = 'Z:\Radar\20140330.007\';
dsearch  = '*.dt0.h5';
files = dir([ddir dsearch]);

rngLims = [60e3 120e3];
Nfft = 256;
Mpul = 1700; % expected number of pulses per beam in a record
dPul = 1700;

ibms = [64982 57281]; Nbms = length(ibms);

fileOpened = 0;
ifile = 1;

while fileOpened == 0    
    fname = files(ifile).name;
    fprintf('File %s\n', fname);
    tm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/RadacTime');
    timeMin = min(min(tm));
    timeMax = max(max(tm));
    disp('Desired Local Time:')
    disp(datestr(utc2date(time,timezone)))
    disp('File Start Local Time:')
    disp(datestr(utc2date(timeMin,timezone)))
    disp('File End Local Time:')
    disp(datestr(utc2date(timeMax,timezone)))
    disp('')
    
    if time > timeMin && time < timeMax
        fileOpened = 1;
        rawtx = hdf5read([ddir fname],'/Raw11/TxData/Samples/Data'); rawtx = squeeze(rawtx(1,:,:,:) + 1.0j*rawtx(2,:,:,:)); rawtx = rawtx(20:end,:,:);
        rawdat = hdf5read([ddir fname],'/Raw11/Data/Samples/Data'); rawdat = squeeze(rawdat(1,:,:,:) + 1.0j*rawdat(2,:,:,:));
        rawbm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/BeamCode');
        rng = hdf5read([ddir fname],'/Raw11/Data/Samples/Range');
        
        Ntx = size(rawtx,1);
        [Nranges, Npulses, Nrecs] = size(rawdat);
        
        Irng = find(rng>rngLims(1) & rng<rngLims(2));
        rng2 = rng(Irng);
    else
        % jump to the correct file
        % find the time increment per file
        deltaTimePerFile = timeMax - timeMin;
        deltaTime = time - timeMin;
        numFilesToJump = floor(deltaTime/deltaTimePerFile);
        ifile = ifile + numFilesToJump;
    end
end

%find the record where the time is
recordOpened = 0;
irec = 1;

while recordOpened == 0 && irec < 10
    timeMin = min(tm(:,irec));
    timeMax = max(tm(:,irec));
    
    if time > timeMin && time < timeMax
        recordOpened = 1;
    else
        irec = irec + 1;
    end
    timeMin = min(tm(:,irec));
    timeMax = max(tm(:,irec));
end

%we have the record and file, now let's plot the data
tbm = rawbm(:,irec);
odat = zeros(Nbms,length(Irng),Mpul);
tms = zeros(Nbms,Mpul);

for ibm = 1:Nbms
    Ibm = find(tbm==ibms(ibm));
    
    tdat = rawdat(:,Ibm,irec);
    ttx = rawtx(:,Ibm,irec);
    tms(ibm,1:length(Ibm)) = tm(Ibm,irec);
    
    Npul = size(tdat,2);
    
    for ipul = 1:Npul
        
        fprintf('\t rec %d/%d, bm %d/%d, pulse %d/%d\n',irec,Nrecs,ibm,Nbms,ipul,Npul);
        
        for irng = 1:length(Irng)
            if Irng(irng)+Ntx-1 <= length(rng)
                x = conj(ttx(:,ipul)).*tdat(Irng(irng):(Irng(irng)+Ntx-1),ipul);
                tpsd = fftshift(abs(fft(x,Nfft)).^2);
                odat(ibm,irng,ipul) = max(tpsd);
            end
        end
    end
end

%%% Plotting

tstr = sprintf('Seconds after %s',datestr(utc2date(tm(1,irec),timezone)));

%for ipul = 1:(dPul/2):(Mpul-dPul)

%ibeg = ipul;
%iend = min([ipul+dPul,Mpul]);
ibeg = 1;
iend = Mpul;

figure

for ibm = 1:Nbms
    
    x = squeeze(odat(ibm,:,ibeg:iend));
    %x = medfilt2(x,[3,3]);
    
    subplot(2,1,ibm);
    imagesc(tms(ibm,ibeg:iend)-tm(1,irec),rng2/1e3,10*log10(x./median(median(x))));
    set(gca,'ydir','normal');
    ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
    colorbar;
end

subplot(2,1,2);
xlabel(tstr);
disp('Desired Local Time:')
disp(datestr(utc2date(time,timezone)))
disp('Record Start Local Time:')
disp(datestr(utc2date(timeMin,timezone)))
disp('Record End Local Time:')
disp(datestr(utc2date(timeMax,timezone)))
disp('')
disp('Click on the coordinates of the head')
headCoordinates = ginput(2);

%get data from head
narrowBeam = squeeze(odat(1,:,ibeg:iend))';
wideBeam = squeeze(odat(2,:,ibeg:iend))';
xStartTime = min(headCoordinates(:,1));
xEndTime = max(headCoordinates(:,1));
yStartHeight = min(headCoordinates(:,2));
yEndHeight = max(headCoordinates(:,2));
[dummy, xStartIndex] = min(abs(tms(ibm,ibeg:iend)-tm(1,irec)-xStartTime));
[dummy, xEndIndex] = min(abs(tms(ibm,ibeg:iend)-tm(1,irec)-xEndTime));
[dummy, yStartIndex] = min(abs(rng2/1e3-yStartHeight));
[dummy, yEndIndex] = min(abs(rng2/1e3-yEndHeight));

%go through each time profile and find the maximum power of the head echo
headPowerNarrowBeam = zeros(xEndIndex-xStartIndex+1,1);
headRangeNarrowBeam = zeros(xEndIndex-xStartIndex+1,1);
headPowerWideBeam = zeros(xEndIndex-xStartIndex+1,1);
headRangeWideBeam = zeros(xEndIndex-xStartIndex+1,1);
for xIndexCurrent = xStartIndex:xEndIndex
    [powerNB, indexNB] = max(narrowBeam(xIndexCurrent,yStartIndex:yEndIndex));
    [powerWB, indexWB] = max(wideBeam(xIndexCurrent,yStartIndex:yEndIndex));
    headPowerNarrowBeam(xIndexCurrent-xStartIndex+1) = powerNB;
    headRangeNarrowBeam(xIndexCurrent-xStartIndex+1) = rng2(indexNB)/1e3;
    headPowerWideBeam(xIndexCurrent-xStartIndex+1) = powerWB;
    headRangeWideBeam(xIndexCurrent-xStartIndex+1) = rng2(indexWB)/1e3;
end

titleFont = 20;
axisFont = 16;
tickFont = 15;

figure
imagesc(tms(1,ibeg+xStartIndex-1:ibeg+xEndIndex-1)-tm(1,irec),rng2/1e3,10*log10(squeeze(odat(1,:,ibeg+xStartIndex-1:ibeg+xEndIndex-1))./median(median(odat(1,:,ibeg+xStartIndex-1:ibeg+xEndIndex-1)))));
set(gca,'ydir','normal');
set(gca,'FontSize',tickFont)
ylabel(sprintf('Range (km)'),'fontsize',axisFont)
xlabel(tstr, 'fontsize',axisFont);
title('Narrow Beam RTI','fontsize',titleFont)
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','SNR (dB)','fontsize',axisFont)
set(cbar_handle,'FontSize',tickFont)

figure
imagesc(tms(2,ibeg+xStartIndex-1:ibeg+xEndIndex-1)-tm(2,irec),rng2/1e3,10*log10(squeeze(odat(2,:,ibeg+xStartIndex-1:ibeg+xEndIndex-1))./median(median(odat(2,:,ibeg+xStartIndex-1:ibeg+xEndIndex-1)))));
set(gca,'ydir','normal');
set(gca,'FontSize',tickFont)
ylabel(sprintf('Range (km)'),'fontsize',axisFont)
xlabel(tstr, 'fontsize',axisFont);
title('Wide Beam RTI','fontsize',titleFont)
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','SNR (dB)','fontsize',axisFont)
set(cbar_handle,'FontSize',tickFont)
    
%end
%save data
save(['C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\Meteors\' meteor_name '.mat'],'ibeg','iend','irec','odat','rng2','tm','tms')
