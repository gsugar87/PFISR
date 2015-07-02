function [ gain ] = radarGainFromAzEl( az,el,beam )
%radarGainFromAzEl For a given azimuth and elevation, this function returns
%the PFISR radar gain for the wide beam (beam = 1) or narrow beam (beam =
%2)
%   Detailed explanation goes here

%set needed variables
c = 2.99792458e8;            %speed of light (m/s)
freq = 450e6;
wavelen = c/freq;
k = 2*pi/wavelen;
% set rotation degrees of the PFISR array
rotz = 14;
roty = 16;
rzMatrix = [[cosd(rotz) -sind(rotz) 0];[sind(rotz) cosd(rotz) 0];[0 0 1]];
ryMatrix = [[cosd(roty) 0 sind(roty)];[0 1 0];[-sind(roty) 0 cosd(roty)]];
rotationMatrix = rzMatrix*ryMatrix;
rxpp = rotationMatrix*[1;0;0];
rypp = rotationMatrix*[0;1;0];
%get the unit vector elements for the el and az we want
rx = cosd(el).*cosd(az);
ry = cosd(el).*sind(az);
rz = sind(el);
%find the x and y angles needed to get the beam gain
xAngle = pi/2-acos(rxpp(1).*rx+rxpp(2).*ry+rxpp(3).*rz);
yAngle = pi/2-acos(rypp(1).*rx+rypp(2).*ry+rypp(3).*rz);

%% READING IN THE RADAR PHASING DATA

disp('Reading in the antenna element coordinate data')
aeu_coordinates = csvread('\Users\Glenn\Documents\MATLAB\PFISR\RadarPhasingInfo\nec_phases_cal11.txt',0,1);
x = aeu_coordinates(:,3);
y = aeu_coordinates(:,4);
num = aeu_coordinates(:,1);
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
%% END OF READING IN THE PHASING DATA
disp('Generating the radar beam pattern for the given az and el.')
F = sum(exp(1j*(k*(x*xAngle + y*yAngle)+wide_phases)));
%gain = F.*conj(F);
gain = F;
end

