
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
% ifileArray = [];
% irecArray = [];
% ibegArray = [];
% iendArray = [];
% snrNarrow5t10Array = [];
% snrWide5t10Array = [];
% snrNarrow10t15Array = [];
% snrWide10t15Array = [];
% snrNarrow15t20Array = [];
% snrWide15t20Array = [];
% snrNarrow20t25Array = [];
% snrWide20t25Array = [];
% snrNarrowPlus25Array = [];
% snrWidePlus25Array = [];

%for ifile = 1:length(files)
for ifile = ifileArray(end):length(files)
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
            xNarrow = squeeze(odat(1,:,ibeg:iend));
            xWide = squeeze(odat(2,:,ibeg:iend));
            snrNarrow = 10*log10(xNarrow./median(median(xNarrow)));
            snrWide = 10*log10(xWide./median(median(xWide)));
            
            ifileArray = [ifileArray,ifile];
            irecArray = [irecArray,irec];
            ibegArray = [ibegArray,ibeg];
            iendArray = [iendArray,iend];
            snrNarrow5t10Array = [snrNarrow5t10Array, sum(sum((snrNarrow >= 5).*(snrNarrow <= 10)))];
            snrWide5t10Array = [snrWide5t10Array, sum(sum((snrWide >= 5).*(snrWide <= 10)))];
            snrNarrow10t15Array = [snrNarrow10t15Array, sum(sum((snrNarrow >= 10).*(snrNarrow <= 15)))];
            snrWide10t15Array = [snrWide10t15Array, sum(sum((snrWide >= 10).*(snrWide <= 15)))];
            snrNarrow15t20Array = [snrNarrow15t20Array, sum(sum((snrNarrow >= 15).*(snrNarrow <= 20)))];
            snrWide15t20Array = [snrWide15t20Array, sum(sum((snrWide >= 15).*(snrWide <= 20)))];
            snrNarrow20t25Array = [snrNarrow20t25Array, sum(sum((snrNarrow >= 20).*(snrNarrow <= 25)))];
            snrWide20t25Array = [snrWide20t25Array, sum(sum((snrWide >= 20).*(snrWide <= 25)))];
            snrNarrowPlus25Array = [snrNarrowPlus25Array, sum(sum(snrNarrow >= 25))];
            snrWidePlus25Array = [snrWidePlus25Array, sum(sum(snrWide >= 25))];
        end
        close all;
        
    end
save('20140331.005.mat','ifileArray','irecArray','ibegArray','iendArray','snrNarrow5t10Array','snrWide5t10Array','snrNarrow10t15Array','snrWide10t15Array','snrNarrow15t20Array','snrWide15t20Array','snrNarrow20t25Array','snrWide20t25Array','snrNarrowPlus25Array','snrWidePlus25Array')    
end
sprintf('Done with all files!')