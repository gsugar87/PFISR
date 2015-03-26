clear all
close all

% SCRIPT TO PLOT METEOR STREAK HAPPENED 
% ON THE 30th of March 2014
% AT 11:00:15,480
clear all

% ULTRA
    % Grab raw data centered at center of event,
    [data_baseline,~,tUT] = rawDMCreader('E:\PFISR Images\UltraPFRR\2014-03-30\2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,42927:43728,0,[100,1100],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,43322:43335,0,[100,1100],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))

% IXON CAMERA
    % Grab raw data centered at center of event,
    
    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,4008:4409,0,[980, 1070],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,4204:4213,0,[980, 1070],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))


%
% SCRIPT TO PLOT METEOR STREAK HAPPENED 
% ON THE 30th of March 2014
% AT 11:24:40,321
clear all

% ULTRA

    % Grab raw data centered at center of event,
    % The grabbed data account for 15sec of signals between frame 195909-196709
    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,120550:121350,0,[100,1100],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,120924:120976,0,[100,1100],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))

    % Grab raw data centered at center of event,

    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,195909:196709,0,[100,1100],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);

% CLASSIC
    
    % Grab raw data centered at center of event,
    
    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,52337:52738,0,[980, 1070],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,52521:52554,0,[980, 1070],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))

%
%
% SCRIPT TO PLOT METEOR STREAK HAPPENED 
% ON THE 30th of March 2014
% AT 11:48:21,811

% ULTRA CAMERA
    % Grab raw data centered at center of event,
    % The grabbed data account for 15sec of signals between frame 195909-196709
    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,195909:196709,0,[100,1100],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    % The grabbed data account for 491ms from frame 196296-196322
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-03-30/2014-03-30T10-46-CamSer7196.DMCdata',512,512,1,1,196296:196322,0,[100,1100],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))

% IXON CAMERA
      % Grab raw data centered at center of event,
    
    [data_baseline,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,99259:99660,0,[980, 1070],'auto','auto');

    % Doing the sum of the baseline to average out noise
    sum_data_baseline=sum(data_baseline,3);

    % Doing average of the baseline to average out noise
    mean_data_baseline=mean(data_baseline,3);

    % Computing the standard deviation of the baseline
    std_data_baseline=std(single(data_baseline),0,3);


    % Grab data regarding only the meteor strike
    
    [data_streak,~,tUT] = rawDMCreader('/media/ExtBook/Ixon/2014-03-30/2014-03-30T10-58-CamSer1387.DMCdata',512,512,1,1,99451:99468,0,[980, 1070],'auto','auto');

    % Doing sum of the meteor data to see the streak
    sum_data_streak=sum(data_streak,3);

    %Computing the mean for the meteor streak
    mean_data_streak=mean(data_streak,3);

    %Computing the standard deviation for the meteor striek
    std_data_streak=std(single(data_streak),0,3);

    %Possible compuations to show streaks

    %Ratio of mean
    figure;imagesc((mean_data_streak./mean_data_baseline))

    %Ratio of standard deviations
    figure;imagesc((std_data_streak./std_data_baseline))
