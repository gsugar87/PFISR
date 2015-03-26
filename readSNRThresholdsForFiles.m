%This script will read a .mat file that contains the number of snr values
%above thresholds for every file and record in the specified .mat file

load('20140331.005.mat')

%get the file and records where there is high snr
goodIndecesNarrow = snrNarrowPlus25Array > 0;
goodIndecesWide = snrWidePlus25Array > 0;

goodIFilesNarrow = ifileArray(goodIndecesNarrow);
goodIRecNarrow = irecArray(goodIndecesNarrow);
goodIBegNarrow = ibegArray(goodIndecesNarrow);
goodIEndNarrow = iendArray(goodIndecesNarrow);
goodIFilesWide = ifileArray(goodIndecesWide);
goodIRecWide = irecArray(goodIndecesWide);
goodIBegWide = ibegArray(goodIndecesWide);
goodIEndWide = iendArray(goodIndecesWide);

%ddir = '/Volumes/AMISR_019/Data AMISR Poker/20140327.017/';
ddir = 'F:\20140331.005\';
%ddir = 'C:\PFISR\Data\';
dsearch  = '*.dt0.h5';
files = dir([ddir dsearch]);

rngLims = [70e3 115e3];
Nfft = 256;
Mpul = 1700; % expected number of pulses per beam in a record
dPul = 300;

ibms = [64982 57281]; Nbms = length(ibms);

load('20140331.005.mat')

for ifileIndex = 9:2:length(goodIFilesNarrow)
    ifile = goodIFilesNarrow(ifileIndex);
    irec = goodIRecNarrow(ifileIndex);
    ibeg = goodIBegNarrow(ifileIndex);
    iend = goodIEndNarrow(ifileIndex);
    
    fname = files(ifile).name;
    fprintf('File %s\n', fname);
    
    rawtx = hdf5read([ddir fname],'/Raw11/TxData/Samples/Data'); rawtx = squeeze(rawtx(1,:,:,:) + 1.0j*rawtx(2,:,:,:)); rawtx = rawtx(20:end,:,:);
    rawdat = hdf5read([ddir fname],'/Raw11/Data/Samples/Data'); rawdat = squeeze(rawdat(1,:,:,:) + 1.0j*rawdat(2,:,:,:));
    rawbm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/BeamCode');
    rng = hdf5read([ddir fname],'/Raw11/Data/Samples/Range');
    tm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/RadacTime');
    
    Ntx = size(rawtx,1);
    [Nranges, Npulses, Nrecs] = size(rawdat);
    
    Irng = find(rng>rngLims(1) & rng<rngLims(2));
    rng2 = rng(Irng);
    figure();
    %read the irec from the file
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
                
                x = conj(ttx(:,ipul)).*tdat(Irng(irng):(Irng(irng)+Ntx-1),ipul);
                tpsd = fftshift(abs(fft(x,Nfft)).^2);
                odat(ibm,irng,ipul) = max(tpsd);
            end
            
        end
    end
    
    %%% Plotting
    tstr = sprintf('Seconds after %s',datestr(utc2date(tm(1,irec),-8)));    %***USED TO BE -9 INSTEAD OF -8 (CHANGED ON 8/8/2014)
    xNarrow = squeeze(odat(1,:,ibeg:iend));
    xWide = squeeze(odat(2,:,ibeg:iend));
    snrNarrow = 10*log10(xNarrow./median(median(xNarrow)));
    snrWide = 10*log10(xWide./median(median(xWide)));
    
    figure; set(gcf)
    
    subplot(2,1,1);
    imagesc(tms(1,ibeg:iend)-tm(1,irec),rng2/1e3,snrNarrow);
    colorbar;
    set(gca,'ydir','normal');
    ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
    subplot(2,1,2);
    imagesc(tms(2,ibeg:iend)-tm(1,irec),rng2/1e3,snrWide);
    set(gca,'ydir','normal');
    ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
    xlabel(tstr);
    colorbar;
    close all
    
end
close all;
sprintf('Done with all files!')