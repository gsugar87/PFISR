%% test read power aeu files
directory = 'C:\Users\Glenn\Documents\MATLAB\PFISR\Radar\aeuPowers\20140330';

files = dir(directory);
files = files(3:end);
numFiles = numel(files);
netPowers = zeros(1,numFiles);
for fileIndex = 1:numFiles
    powerStruct = xml2struct([directory '\' files(fileIndex).name]);
    numPanels = numel(powerStruct.results.panel);
    netPower = 0;
    summedAEUs = 0;
    for i = 1:numPanels
        numAEUs = numel(powerStruct.results.panel{i}.aeu);
        for j = 1:numAEUs
            netPower = netPower + str2num(powerStruct.results.panel{i}.aeu{j}.Attributes.pwatts);
            summedAEUs = summedAEUs+1;
        end
    end 
    netPowers(fileIndex) = netPower
end