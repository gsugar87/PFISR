%                       rcs2q_sph.m
%
% Converts a modeled q to radar RCS using the 3D full-wave solution.
%
% INPUT
% alt = vector of altitudes (km)
% rcs = rcs (dBsm)
% rr = range rate or 3D velocity (km/s);
%
% OUTPUTS
% q = electron line density (m^-1)
%
% S. Close August 2011
% Edited by G. Sugar June 11 to work with PFISR data

function [qvec,nvec] = rcs2q(alt,rcs,rr)
warning off
global CONST EXPWAVEFORM
freq    = 4.493e8;                     %frequency in Hz
freqmhz = freq/1e6;                %frequency in MHz
lambda  = 2.99e8./freq;                         %lambdalength in m
k       = 2*pi./lambda;
omega   = 2*pi*freq;

%Get background electron density
electrons_height = csvread('C:\Users\Glenn\Documents\CEDAR\Electron Density.csv');
elalt   = electrons_height(:,1);    %altitude in km
nb      = electrons_height(:,2)*1e-6;    %electron number density in cm^-3
nbint   = interp1(elalt,nb,alt,'cubic');    %electron number density in cm^-3

%Get the meanfreepath from the internet and use the number density in there
load meanfreepath.txt;    
mfpalt   = meanfreepath(:,1)/1e3;               %altitude...put it into km
mfp      = meanfreepath(:,8);                   %mean free path in m
numdens  = meanfreepath(:,5);                   %number density 
colfreq  = meanfreepath(:,7);                  %collision frequency

nd       = interp1(mfpalt,numdens,alt,'cubic'); %in m^-3
mfprad   = interp1(mfpalt,mfp,alt,'cubic');
nu       = interp1(mfpalt,colfreq,alt,'cubic');

%Gaussian falls off twice as fast, so where eps = 0 is twice as far out so need half the radius of the 
%Herlofson (parabolic) approximation (if parabolic) --> smaller the radius, smaller the density
inirad   = (.023)*(2.845e18.*(abs(rr)).^(.8))./nd; %(.023 is about 1mfp)  %.022 works well!
%inirad      = (abs(rr).^(1/2)./nu.^(1/2));

%GET REFLECTION COEFFICIENT AS A FUNCTION OF RADIUS AND DENSITY IN 3D
%gamma   = 2*pi./inirad;   
nm      = 10.^(9:.1:21);        %electron density (m^-3)
n       = nm./1e6;              %electron density (cm^-3)
omegac  = 2*pi*9000*sqrt(n);    %plasma frequency

%SET UP LOOP
pl      = [];  ref = [];  sig = [];
%fprintf('Max radius is %3.2fm (wave=%3.2f m).\n',max(inirad),lambda);
for ii = 1:length(inirad)
    %s_top   = log10(log(nmeter./nbint(ii)));
    %WHEN DOES THIS METHOD BREAK DOWN?  WHEN FAILS > 1
    %dr          = inirad(ii)./100;
    %rins        = (0:dr:inirad(ii));
    %rmatins     = (rins'*ones(1,length(omegac)))';
    %omegacins   = omegac'*ones(1,length(rins));
    %eins        = 1 - (omegacins.^2./omega.^2).*exp(-(rmatins./inirad(ii)).^2);
    %fails       = abs(eins.*(k.*rmatins).^2);    %where this is greater than m, my theory fails
    m           = 0;
    gvec        = zeros(1,length(n));
    diffla      = 10*ones(1,length(n));
    while max(abs(diffla))>2 && m < 4
    %for m = 0:0
        m               = m+1;
        %figure(1);imagesc(rins,nm,fails,[0 m*(m+1)]);axis xy;colorbar;xlabel('R out to head echo radius');ylabel('Plasma frequency');
        %[rovec,A,B]     = headiterate_varyingrstop(omega,n,inirad(ii),nbint(ii),m);
        [rovec,A,B]     = headiterate(omega,n,inirad(ii),nbint(ii),m);
        %r               = rovec.*inirad(ii);
        r               = inirad(ii);
        %r               = 1;
        J               = besselj(m,k.*r);
        H2              = besselh(m,2,k.*r);
        %aminv           = -(2-(m.*H2.*A.*(r.^(2*m+1))./(B.*J.*(m+1))));    %MINE
        aminv           = -(1+((i.*m.*factorial(2.*m).*factorial(2.*m+1).*A)./((4.^m).*(m+1).*(factorial(m).^2).*(k.^(2*m+1)).*B)));
        am              = abs(1./aminv);
        %gvecam          = ((m+(1/2)).*am);  
        %sigmaam         = db((gvecam.^2.*lambda.^2./pi));
        %figure(12);semilogx(nm,sigmaam);keyboard
        last            = gvec;
        gvec            = ((m+(1/2)).*am) + last;
        diffla          = last - gvec;  %difference between last gvec and this one
        %fprintf('ii = %d, m = %d, max(diffla) = %3.2d\n',ii,m,max(abs(diffla)));
    end
    %fprintf('Summed m to %3.2d,\n',m);
    sigma   = db((gvec.^2.*lambda.^2./pi));
    diffv   = abs(rcs(ii)-sigma);
    [ix,ixcc]   = min(diffv);
    md = find(diffv < 1);   %takes care of Mie resonance, finds all differences < 1
    if length(md) > 1, ixcc = md(1); end;
    ref = [ref; am(ixcc)];  
    sig = [sig; sigma(ixcc)];
    pl = [pl; nm(ixcc)];
    gvec = 0;
    %semilogx(nm,sigma);xlabel('Modeled Plasma Density (el/m^3)');ylabel('Modeled RCS (dBsm)');grid on
    if 1 == 0
        if freq == 422
            figure(2);semilogx(nm,sigma,'c',nm,rcs(ii),'b');hold off;
        else
            rcsvec = rcs(ii).*ones(1,length(nm));
            figure(2);h=semilogx(nm,sigma,'r',nm,rcsvec,'k');
            set(h(1),'linewidth',3);
        end
        keyboard
    end
    %keyboard
end
nvec = pl;
minnvec = nvec.*exp(-(1)^2);       %evaluate at r = rmax
qvec = minnvec.*pi.*(inirad.^2);
%keyboard

%global nvec inirad
%num_electrons = quad(@nv,0,inirad,1e-6);
%qvec = num_electrons./(2.*inirad);

%rod = (bd(rcs).*(lambda.^4)./(64*(pi^5))).^(1/6);  %get radius from RCS using overdense in m
%qod = 3.16e14.*bd(rcs);
%semilogy(alt,nvec,'*',alt,qvec,'*',alt,qod,'g.');legend('n','q','qod');keyboard
%[a,b]=max(nvec);
%fprintf('Max nvec is %3.2d e/m^3 at %3.2f km.\n\n',nvec(b),alt(b));
%keyboard
%figure(4);semilogy(2*pi*inirad./lambda,nvec,'.');
%figure(4);semilogy(alt,nvec,'g.');