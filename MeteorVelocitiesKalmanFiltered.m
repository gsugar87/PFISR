function [trackedLocation] = MeteorVelocitiesKalmanFiltered(input_positions)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
trackedLocation=input_positions;
 kalmanFilter = []; isTrackInitialized = false;
 for idx = 1: size(input_positions,1);
     detectedLocation = input_positions(idx,:);
     isObjectDetected = ~isnan(detectedLocation(1,1));

     if ~isTrackInitialized
       if isObjectDetected
         kalmanFilter = configureKalmanFilter('ConstantAcceleration',detectedLocation(1,:), [1 1 1]*1e3, [15, 15, 15], 15);
         isTrackInitialized = true;
       end
     else
       if isObjectDetected
         predict(kalmanFilter);
         trackedLocation(idx,:) = correct(kalmanFilter, detectedLocation(1,:));
       else
         trackedLocation(idx,:) = predict(kalmanFilter);
       end
     end
 end

end