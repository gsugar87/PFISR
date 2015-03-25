function time_ml_format=utc2date(UTCtime,Timezone)
% converts the unix time to matlab time

seconds_per_day = 3600*24;
seconds_per_hour = 60*60; 

% Matlab has 1-Jan-0000 as reference day and utc is in seconds,
% here we calculate time from 0000 to 1904

% Calculating the matlab format
time_ml_format=(UTCtime + Timezone*seconds_per_hour)/(seconds_per_day) + datenum(1970,1,1);