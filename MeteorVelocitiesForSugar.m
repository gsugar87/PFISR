close all
clear all

tic
%filename = 'C:/Users/Lore/Documents/PhD/PokerFlatMar2014/31-03-2014DetectedEvents.xlsx';
filename = 'C:/Users/Glenn/Documents/MATLAB/PFISR/Spreadsheets/31-03-2014DetectedEvents20150713.xlsx';
%
% Call MeteorVelocitiesTracking
% Call MeteorVelocitiesKalmanFiltered
% Call MeteorVelocitiesIntensity
% Call FindPixelsOfInterestForMeteor
% Call Bresenham
%
NightToBeAnalyzed='F:/2014-03-31/2014-03-31T06-12-CamSer7196.DMCdata';
%TrueMeteorBlobs='C:\Users\Lore\Documents\PhD\PokerFlatMar2014\31032014TrueMeteorsBlobs.xlsx';
TrueMeteorBlobs='C:/Users/Glenn/Documents/MATLAB/PFISR/Spreadsheets/31032014TrueMeteorsBlobs.xlsx';
MeteorEvents = xlsread(filename);
location=[];
fittedlocation=[];
kalmanlocation=[];
dsearch  = '*.mat';
dirofthing = 'E:\SavedEvents31st\';
file=dir(['E:\SavedEvents31st\' dsearch])
for kkk=1:size(file,1)
    if kkk>=74
        kkk=kkk+1;
        load([dirofthing file(kkk-1).name])
    else
        load([dirofthing file(kkk).name])
    end
    if MeteorEvents(kkk,11)==1
            fprintf('Getting camera pixel location\n');
    [output_args] = MeteorVelocitiesTracking(NightToBeAnalyzed,data_event,...
        MeteorEvents(kkk,1),MeteorEvents(kkk,2),MeteorEvents(kkk,6),MeteorEvents(kkk,7),...
        MeteorEvents(kkk,8),MeteorEvents(kkk,9),TrueMeteorBlobs,data_baseline);
    % kalman filtering for filling the blanks if object goes missing
    [kalman_filtered_locations] = MeteorVelocitiesKalmanFiltered(output_args(:,1:2));
    X=kalman_filtered_locations(:,1);
    [P,S] = polyfit(X,kalman_filtered_locations(:,2),2);
    [Y,DELTA] = polyval(P,kalman_filtered_locations(:,1),S);
    % [P1,S1] = polyfit(kalman_filtered_locations(:,2),kalman_filtered_locations(:,1),2);
    % [X,DELTA1] = polyval(P1,kalman_filtered_locations(:,2),S1);
    fprintf('Getting Meteor Signal data\n');
    [Meteor_Intensity] = MeteorVelocitiesIntensity( kalman_filtered_locations,data_event,data_baseline, var_baseline );
    location=[location;output_args];
    kalmanlocation=[kalmanlocation;kalman_filtered_locations];
    fittedlocation=[fittedlocation;X,Y];
    figure(100)
    imagesc(data_event(:,:,1))
    hold on
    plot(kalmanlocation(:,1),kalmanlocation(:,2),'ko',location(:,1),location(:,2),'ro', fittedlocation(:,1),fittedlocation(:,2),'go')
%   
    pause
    close all
    end


% end

% pause

end