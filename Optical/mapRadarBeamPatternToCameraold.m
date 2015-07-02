% This script will read in the radar beam parameters and generate the beam
% pattern
%this uses the model given by equation 1 in Chau et al 2009
% F(thetax,thetay) = sum_overi(gi*exp(j*k*(xi*thetax+yi*thetay) + j*phii))

%set which format we want to generate the beam pattern
%1: use linspace thetax thetay (no mapping to camera)
%2: use camera pixels to get thetax and thetay (mapping to camera)
mapType = 2;
%mapType = 2;

%set needed variables
freq = 450e6;
wavelen = 3e8/freq;
numPixels = 512;
centerPixelCoords = [271,246];
degreesPerPixel = 0.0177;
% set rotation degrees of the PFISR array
rotz = 14;
roty = 16;
rzMatrix = [[cosd(rotz) -sind(rotz) 0];[sind(rotz) cosd(rotz) 0];[0 0 1]];
ryMatrix = [[cosd(roty) 0 sind(roty)];[0 1 0];[-sind(roty) 0 cosd(roty)]];
rotationMatrix = rzMatrix*ryMatrix;
rxpp = rotationMatrix*[1;0;0];
rypp = rotationMatrix*[0;1;0];

if mapType == 2
    %create theta arrays to map to camera coordinates
    %use the calibration data from the cameras
    disp('loading the camera azimuth and elevation data')
    load('ixonCalibrate2.mat')
    azC = az;
    azC(find(azC > 200)) = azC(find(azC > 200)) - 360;
    elC = el;
    load('ultraCalibrate2.mat')
    azU = az;
    azU(find(azU > 200)) = azU(find(azU > 200)) - 360;
    elU = el;
    elCenter = 74;
    azCenter = 14;
    classicAnglesX = cos(elC*pi/180).*(azC-azCenter);
    classicAnglesY = elC-elCenter;
    ultraAnglesX = cos(elU*pi/180).*(azU-azCenter);
    ultraAnglesY = elU-elCenter;
    % phi = atan2((elU-elCenter).*pi./180.*sin(elCenter.*pi./180),(azU-azCenter).*cos(elCenter.*pi./180));
    % uaX = (elU-elCenter).*sin(phi);
    % uaY = (elU-elCenter).*cos(phi);
    % uaX = (azU-azCenter).*cos(elU*pi/180);
    % uaY = elU-elCenter;
    %convert to radians
    ultraAnglesX = ultraAnglesX*pi/180;
    ultraAnglesY = ultraAnglesY*pi/180;
    classicAnglesX = classicAnglesX*pi/180;
    classicAnglesY = classicAnglesY*pi/180;

    %% TESTING%%
    %lat = elCenter+90;
    %convert the azimuth and elevation to radar beam az el using 3 formulas
    %below:
    %dec=asin(sin(el)*sin(lat)+cos(el)*cos(lat)*cos(az))
    %LHA = atan2(-sin(az)*cos(el)/cos(dec),
    %(sin(el)-sin(dec)*sin(lat))/(cos(dec)*cos(lat)))
    %RA = mod(theta_localSolarTime-LHA,360) (this is not used)
    disp('converting azimuth and elevation to radar beam az el')
    lat = elCenter-90;
    decU = asind(sind(elU).*sind(lat)+cosd(elU).*cosd(lat).*cosd(azU-azCenter));
    LHAU = atan2(-sind(azU-azCenter).*cosd(elU)./cosd(decU), ...
        (sind(elU)-sind(decU).*sind(lat))./(cosd(decU).*cosd(lat))).*(180/pi);
    decC = asind(sind(elC).*sind(lat)+cosd(elC).*cosd(lat).*cosd(azC-azCenter));
    LHAC = atan2(-sind(azC-azCenter).*cosd(elC)./cosd(decC), ...
        (sind(elC)-sind(decC).*sind(lat))./(cosd(decC).*cosd(lat))).*(180/pi);
    txU = decU.*sind(LHAU);
    tyU = decU.*cosd(LHAU);
    txC = decC.*sind(LHAC);
    tyC = decC.*cosd(LHAC);
    %%
    % ultraAnglesX = decU*pi/180;
    % ultraAnglesY = LHAU*pi/180;
    % classicAnglesX = decC*pi/180;
    % classicAnglesY = LHAC*pi/180;
    ultraAnglesX = txU*pi/180;
    ultraAnglesY = tyU*pi/180;
    classicAnglesX = txC*pi/180;
    classicAnglesY = tyC*pi/180;
    %try to get thetax, thetay on 2/3/2015
    ultraAnglesX = pi/2-acos(sind(decU));
    ultraAnglesY = pi/2-acos(cosd(decU).*sind(LHAU));
    %get thetax,theta y on 2/4/2015
    rx = sind(decU).*cosd(LHAU);
    ry = sind(decU).*sind(-LHAU);
    rz = cosd(decU);
    ultraAnglesX = pi/2 - acos(rx.*rotationMatrix(1,1) + ry.*rotationMatrix(1,2) + rz.*rotationMatrix(1,3));
    ultraAnglesY = pi/2 - acos(rx.*rotationMatrix(2,1) + ry.*rotationMatrix(2,2) + rz.*rotationMatrix(2,3));
    %get thetax,theta y on 2/4/2015 (again)
    rx = cosd(elU).*cosd(azU);
    ry = cosd(elU).*sind(azU);
    rz = sind(elU);
    ultraAnglesX = pi/2-acos(rxpp(1).*rx+rxpp(2).*ry+rxpp(3).*rz);
    ultraAnglesY = pi/2-acos(rypp(1).*rx+rypp(2).*ry+rypp(3).*rz);
    
    disp('Generating the thetax thetay arrays')
    thetax = linspace(1-centerPixelCoords(1),numPixels-centerPixelCoords(1),numPixels);
    thetay = linspace(1-centerPixelCoords(2),numPixels-centerPixelCoords(2),numPixels);
    thetax = thetax*pi/180*degreesPerPixel;
    thetay = thetay*pi/180*degreesPerPixel;

    thetax = linspace(1-centerPixelCoords(1),numPixels-centerPixelCoords(1),numPixels);
    thetay = linspace(1-centerPixelCoords(2),numPixels-centerPixelCoords(2),numPixels);
    thetax = thetax*pi/180*degreesPerPixel;
    thetay = thetay*pi/180*degreesPerPixel;
elseif mapType == 1
    thetax = linspace(-4.5,4.5,512);
    thetay = linspace(-4.5,4.5,512);
    thetax = thetax*pi/180;
    thetay = thetay*pi/180;
end

k = 2*pi/wavelen;
disp('Reading in the antenna element coordinate data')
aeu_coordinates = csvread('\Users\Glenn\Documents\MATLAB\PFISR\RadarPhasingInfo\nec_phases_cal11.txt',0,1);
x = aeu_coordinates(:,3);
y = aeu_coordinates(:,4);
num = aeu_coordinates(:,1);
meanx = mean(x);
meany = mean(y);
fid = fopen('\Users\Glenn\Documents\MATLAB\PFISR\RadarPhasingInfo\nec_phases_cal11.txt');
coord_file = fread(fid,'*char');
fclose(fid);
C_rows = regexp(transpose(coord_file),'R[0-9][0-9]','match');
C_cols = regexp(transpose(coord_file),'C[0-9][0-9]','match');
%loop through each element and make a matrix holding the row, column, number, and
%coordinates (RCXY = row, column, num, x, y)
RCNXY = zeros(4096,5);
for i = 1:4096
    RCNXY(i,:) = [str2num(C_rows{i}(1,2:3)),str2num(C_cols{i}(1,2:3)), ...
        num(i), x(i),y(i)];
end
phix = 0.033;
phiy = 0.029;
para_phases = phix*(x-meanx).^2 + phiy*(y-meany).^2; 
disp('Reading in the antenna element phase data')
fid = fopen('\Users\Glenn\Documents\MATLAB\PFISR\RadarPhasingInfo\phases24513.txt');
W_phases = fread(fid,'*char');
fclose(fid);
fid = fopen('\Users\Glenn\Documents\MATLAB\PFISR\RadarPhasingInfo\phases32214.txt');
N_phases = fread(fid,'*char');
fclose(fid);
%get the rows and columns for the file
N_rows = regexp(transpose(N_phases),'R[0-9][0-9]','match');
N_cols = regexp(transpose(N_phases),'C[0-9][0-9]','match');
N_nums = regexp(transpose(N_phases),'AEU[0-9][0-9]','match');
N_start = regexp(transpose(N_phases),'[^\s]AEU*[0-9][^\s]*[0-9]')+7;
N_end = regexp(transpose(N_phases),'[\n]');
W_rows = regexp(transpose(W_phases),'R[0-9][0-9]','match');
W_cols = regexp(transpose(W_phases),'C[0-9][0-9]','match');
W_nums = regexp(transpose(W_phases),'AEU[0-9][0-9]','match');
W_start = regexp(transpose(W_phases),'[^\s]AEU*[0-9][^\s]*[0-9]')+7;
W_end = regexp(transpose(W_phases),'[\n]');
wide_phases = zeros(4096,1);
narrow_phases = zeros(4096,1);

% %calculate generic thetax and thetay
% angleLim = 10;
% dtheta = 0.001;
% thetax = -angleLim*pi/180:dtheta:angleLim*pi/180+dtheta;
% thetay = thetax;


%loop through each element and make a matrix holding the row, column, and
%phases (RCNP = row, column, number, phase)
narrowRCNP = zeros(4096,4);
wideRCNP = zeros(4096,4);
for i = 1:4096
    narrowRCNP(i,:) = [str2num(N_rows{i}(1,2:3)),str2num(N_cols{i}(1,2:3)), ...
        str2num(N_nums{i}(1,4:5)),str2num(transpose(N_phases(N_start(i):N_end(i))))*pi/180];
    wideRCNP(i,:) = [str2num(W_rows{i}(1,2:3)),str2num(W_cols{i}(1,2:3)), ...
        str2num(W_nums{i}(1,4:5)),str2num(transpose(W_phases(W_start(i):W_end(i))))*pi/180];
end
%loop through each antenna element (again) and find where the phases match
%the coordinates
for i = 1:4096
    narrow_index = find(narrowRCNP(:,1) == RCNXY(i,1) & narrowRCNP(:,2) == RCNXY(i,2) & narrowRCNP(:,3) == RCNXY(i,3));
    wide_index = find(wideRCNP(:,1) == RCNXY(i,1) & wideRCNP(:,2) == RCNXY(i,2) & wideRCNP(:,3) == RCNXY(i,3));
    narrow_phases(i) = narrowRCNP(narrow_index,4);
    wide_phases(i) = wideRCNP(wide_index,4);
end

%Calculate the gains for each beam pattern
disp('Generating the radar beam pattern gains for each camera pixel')
F_Narrow = zeros(numel(thetax),numel(thetay));
F_Wide = zeros(numel(thetax),numel(thetay));
F_Para = zeros(numel(thetax),numel(thetay));
F_N = F_Narrow;
g = ones(1,4096);
for a = 1:numel(thetax)
    %F_Narrow(a,:) = sum(exp(j*(k*(repmat(x*thetax(a),1,length(thetax))+mtimesx(y,thetay))+repmat(narrow_phases,1,length(thetay)))),1);
    %F_Wide(a,:) = sum(exp(j*(k*(repmat(x*thetax(a),1,length(thetax))  +mtimesx(y,thetay))+repmat(wide_phases,1,length(thetay)))),1);
    %F_Para(a,:) = sum(exp(j*(k*(repmat(x*thetax(a),1,length(thetax)) + mtimesx(y,thetay))+repmat(para_phases,1,length(thetay)))),1);
    F_Narrow(a,:) = sum(exp(1j*(k*(mtimesx(x,ultraAnglesX(a,:))+mtimesx(y,ultraAnglesY(a,:)))+repmat(narrow_phases,1,size(ultraAnglesX,1)))),1);
    F_Wide(a,:) = sum(exp(1j*(k*(mtimesx(x,ultraAnglesX(a,:))+  mtimesx(y,ultraAnglesY(a,:)))+repmat(wide_phases,1,size(ultraAnglesX,1)))),1);
    F_Para(a,:) = sum(exp(1j*(k*(mtimesx(x,ultraAnglesX(a,:))+  mtimesx(y,ultraAnglesY(a,:)))+repmat(para_phases,1,size(ultraAnglesX,1)))),1);
%     for b = 1:numel(thetay)
%          %F_N(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+narrow_phases));
%          %F_W(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+wide_phases));
%          %F_P(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+para_phases));
%          F_N(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+narrow_phases));
%          F_W(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+wide_phases));
%          %F_P(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+para_phases));
%     end
%         F_Narrow(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+narrow_phases));
%         F_Wide(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+wide_phases));
%         F_Para(a,b) = g*exp(j*(k*(x*ultraAnglesX(a,b)+y*ultraAnglesY(a,b))+para_phases));
%         [F_Narrow(a,b),F_Wide(a,b)] = radarGain(thetax(a),thetay(b));
    %end
end

gainNarrow = db(F_Narrow.*conj(F_Narrow));
gainWide = db(F_Wide.*conj(F_Wide));
gainPara = db(F_Para.*conj(F_Para));
gainNarrow = gainNarrow - max(max(gainNarrow));
gainWide = gainWide - max(max(gainWide));
gainPara = gainPara - max(max(gainPara));

%Plotting...
thetaxDeg = thetax*180/pi;
thetayDeg = thetay*180/pi;
[C,h] = contour(thetaxDeg,thetayDeg,gainNarrow,[-40 -30 -20 -10 -3]);
hcl = clabel(C,h,'FontSize',10,'Color','k','Rotation',0);
title('Narrow Beam')
xlabel('theta x')
ylabel('theta y')

[C,h] = contour(thetaxDeg,thetayDeg,gainWide,[-40 -30 -20 -10 -3]);
hcl = clabel(C,h,'FontSize',10,'Color','k','Rotation',0);
title('Wide Beam')
xlabel('theta x')
ylabel('theta y')

[C,h] = contourf(thetaxDeg,thetayDeg,gainWide,[-40 -30 -20 -10 -3]);
colorbar
title('Wide Beam Gain (dB)')
xlabel('\theta_x (degrees)')
ylabel('\theta_y (degrees)')

[C,h] = contourf(thetaxDeg,thetayDeg,gainNarrow,[-40 -30 -20 -10 -3]);
colorbar
title('Narrow Beam Gain (dB)')
xlabel('\theta_x (degrees)')
ylabel('\theta_y (degrees)')

[C,h] = contourf(gainNarrow,[-40 -30 -20 -10 -3]);
colorbar
title('Narrow Beam')
xlabel('X Pixel')
ylabel('Y Pixel')

[C,h] = contourf(gainWide,[-40 -30 -20 -10 -3]);
colorbar
title('Wide Beam')
xlabel('X Pixel')
ylabel('Y Pixel')