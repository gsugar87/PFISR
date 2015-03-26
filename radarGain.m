function [ narrowGain, wideGain ] = radarGain( thetaX, thetaY )
%radarGain This function accepts a thetax and thetay and returns the PFISR
%radar gains at the specified angle
%   Detailed explanation goes here

freq = 450e6;
wavelen = 3e8/freq;
k = 2*pi/wavelen;
aeu_coordinates = csvread('RadarPhasingInfo\nec_phases_cal11.txt',0,2);
x = aeu_coordinates(:,2);
y = aeu_coordinates(:,3);
fid = fopen('RadarPhasingInfo\phases24513.txt');
W = fread(fid,'*char');
fclose(fid);
fid = fopen('RadarPhasingInfo\phases32214.txt');
N = fread(fid,'*char');
fclose(fid);
W_start = regexp(transpose(W),'[^\s]AEU*[0-9][^\s]*[0-9]')+7;
W_end = regexp(transpose(W),'[\n]');
N_start = regexp(transpose(W),'[^\s]AEU*[0-9][^\s]*[0-9]')+7;
N_end = regexp(transpose(W),'[\n]');
wide_phases = zeros(4096,1);
narrow_phases = zeros(4096,1);
for i = 1:4096
    wide_phases(i) = str2num(transpose(W(W_start(i):W_end(i))));
    narrow_phases(i) = str2num(transpose(N(N_start(i):N_end(i))));
end
g = ones(4096,1);

narrowGain = sum(g.*exp(j*(k*(thetaX*x+thetaY*y)+narrow_phases)));
wideGain = sum(g.*exp(j*(k*(thetaX*x+thetaY*y)+wide_phases)));

end

