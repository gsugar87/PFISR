% This script will read in the radar beam parameters and generate the beam
% pattern
%this uses the model given by equation 1 in Chau et al 2009
% F(thetax,thetay) = sum_overi(gi*exp(j*k*(xi*thetax+yi*thetay) + j*phii))

freq = 450e6;
wavelen = 3e8/freq;
numPixels = 512;
centerPixelCoords = [271,246];
degreesPerPixel = 0.0177;

%create theta arrays to map to camera coordinates
thetax = linspace(1-centerPixelCoords(1),numPixels-centerPixelCoords(1),numPixels);
thetay = linspace(1-centerPixelCoords(2),numPixels-centerPixelCoords(2),numPixels);
thetax = thetax*pi/180*degreesPerPixel;
thetay = thetay*pi/180*degreesPerPixel;

k = 2*pi/wavelen;
aeu_coordinates = csvread('RadarPhasingInfo\nec_phases_cal11.txt',0,1);
x = aeu_coordinates(:,3);
y = aeu_coordinates(:,4);
num = aeu_coordinates(:,1);
meanx = mean(x);
meany = mean(y);
fid = fopen('RadarPhasingInfo\nec_phases_cal11.txt');
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
fid = fopen('RadarPhasingInfo\phases24513.txt');
W_phases = fread(fid,'*char');
fclose(fid);
fid = fopen('RadarPhasingInfo\phases32214.txt');
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
F_Narrow = zeros(numel(thetax),numel(thetay));
F_Wide = zeros(numel(thetax),numel(thetay));
F_Para = zeros(numel(thetax),numel(thetay));
g = ones(1,4096);
for a = 1:numel(thetax)
    for b = 1:numel(thetay)
        F_Narrow(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+narrow_phases));
        F_Wide(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+wide_phases));
        F_Para(a,b) = g*exp(j*(k*(x*thetax(a)+y*thetay(b))+para_phases));
%         [F_Narrow(a,b),F_Wide(a,b)] = radarGain(thetax(a),thetay(b));
    end
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
title('Wide Beam')
xlabel('theta x')
ylabel('theta y')

[C,h] = contourf(thetaxDeg,thetayDeg,gainNarrow,[-40 -30 -20 -10 -3]);
colorbar
title('Narrow Beam')
xlabel('theta x')
ylabel('theta y')