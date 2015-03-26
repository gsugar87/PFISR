
%ddir = '/Volumes/AMISR_019/Data AMISR Poker/20140327.017/';
ddir = 'E:\20140331.005\';
%ddir = 'C:\PFISR\Data\';
dsearch  = '*.dt0.h5';
files = dir([ddir dsearch]);

rngLims = [70e3 115e3];
Nfft = 256;
Mpul = 1700; % expected number of pulses per beam in a record
dPul = 300;

ibms = [64982 57281]; Nbms = length(ibms);

for ifile = 1:length(files)
        
    fname = files(ifile).name;
    fprintf('File %s\n', fname);

    %dat = read_data([ddir fname]);
    % 
    %rawtx = dat.Raw11.TxData.Samples.Data; rawtx = squeeze(rawtx(1,:,:,:) + 1.0j*rawtx(2,:,:,:)); rawtx = rawtx(20:end,:,:);
    %rawdat = dat.Raw11.Data.Samples.Data; rawdat = squeeze(rawdat(1,:,:,:) + 1.0j*rawdat(2,:,:,:));
    %rawbm = dat.Raw11.TxData.RadacHeader.BeamCode;
    %rng = dat.Raw11.Data.Samples.Range;
    %tm = dat.Raw11.TxData.RadacHeader.RadacTime;
    
    rawtx = hdf5read([ddir fname],'/Raw11/TxData/Samples/Data'); rawtx = squeeze(rawtx(1,:,:,:) + 1.0j*rawtx(2,:,:,:)); rawtx = rawtx(20:end,:,:);
    rawdat = hdf5read([ddir fname],'/Raw11/Data/Samples/Data'); rawdat = squeeze(rawdat(1,:,:,:) + 1.0j*rawdat(2,:,:,:));
    rawbm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/BeamCode');
    rng = hdf5read([ddir fname],'/Raw11/Data/Samples/Range');
    tm = hdf5read([ddir fname],'/Raw11/TxData/RadacHeader/RadacTime');
    
    Ntx = size(rawtx,1);
    [Nranges, Npulses, Nrecs] = size(rawdat);
    
    Irng = find(rng>rngLims(1) & rng<rngLims(2));
    rng2 = rng(Irng);
    
    for irec = 1:Nrecs
        
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
        
        tstr = sprintf('Seconds after %s',datestr(utc2date(tm(1,irec),-9)));
        
        for ipul = 1:(dPul/2):(Mpul-dPul)
            
            ibeg = ipul;
            iend = min([ipul+dPul,Mpul]);
            
            figure; set(gcf, 'visible', 'off')
            
            for ibm = 1:Nbms
                
                x = squeeze(odat(ibm,:,ibeg:iend));
                %x = medfilt2(x,[3,3]);
                
                subplot(2,1,ibm);
                imagesc(tms(ibm,ibeg:iend)-tm(1,irec),rng2/1e3,10*log10(x));
                set(gca,'ydir','normal');
                ylabel(sprintf('Range (km); bm %d',ibms(ibm)));
                colorbar;
                set(gcf, 'visible', 'off')
            end
            
            subplot(2,1,2);
            xlabel(tstr);
            set(gcf, 'visible', 'off')
            
            oname = sprintf('%s-rec%.3d-%.3d.png',fname,irec,ipul);
            print(gcf,'-dpng','-r150',oname);
            
        end
        close all;
        
    end
    
end