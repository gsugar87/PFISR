function [ output_args ] = MeteorVelocitiesTracking( NightToBeAnalyzed,data_event,EventStart,EventEnd,x1,y1,x2,y2,TrueMeteorBlobs,mean_data_baseline )
%
%
% calls for FindPixelsOfInterestForMeteor
% 
%   Detailed explanation goes here
Location=[];
%Creation of the logit function for blob detection
num = xlsread(TrueMeteorBlobs);
data=[];
for iii=1:size(num,1)
    if isnan(num(iii,1))~=1
    data=[data;num(iii,:)];
    end
end
[b,dev,stats] = glmfit([data(:,1:3),data(:,5)],data(:,4),'binomial','link','logit');
% Finding Blobs
% [data_baseline,~,tUT] = rawDMCreader(NightToBeAnalyzed,512,512,1,1,EventStart-450:EventStart,0,[100,1100],'auto','auto');
%     mean_data_baseline=dataUltraBackground;
% mean_data_baseline=mean(data_baseline,3);
% std_data_baseline=std(double(data_baseline),[],3);
scatto=prctile(mean_data_baseline(:),99); % tells me when probably the pixel is generally bright in the image such that I always discard them when I do tracking
mean_data_bright_pixels=find(mean_data_baseline>scatto);
% Let's look only at the pixels where the meteor goes through
[non_zero_elements_meteora,m,q]=FindPixelsOfInterestForMeteor(x1,y1,x2,y2);
temp1=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
temp2=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
temp1(non_zero_elements_meteora)=mean_data_baseline(non_zero_elements_meteora);
% temp2(non_zero_elements_meteora)=std_data_baseline(non_zero_elements_meteora);
mean_data_baseline=temp1;
% Find bright object 
for zzz=1:size(data_event,3)
    temp1=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
    temp2=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
    temp1=data_event(:,:,zzz);
    temp2(non_zero_elements_meteora)=temp1(non_zero_elements_meteora);
    data_event(:,:,zzz)=temp2;
end

temp1=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
temp2=zeros(512,512); %dummy variable to store non zero mean baseline of interest for the meteor
temp3=temp2;
temp4=temp2;
data_event=double(data_event);

for zzz=1:size(data_event,3)
    temp1=data_event(:,:,zzz)-mean_data_baseline;
%   Since bright pixels have high std and we know where we are, we zero
%   them
    temp1(mean_data_bright_pixels)=0;
    temp1(temp1<0)=0;
    temp1_zero_pixels=find(temp1==0);
    %
    % Apply median and average filter
    %
    temp2=medfilt2(temp1,[3 3]);
    temp2_zero_pixels=find(temp2==0);
    % Apply average filter
    %
    h = fspecial('average',[5 5]);
    temp3= filter2(h,temp2);
    temp3(temp1_zero_pixels)=0;
    temp3(temp2_zero_pixels)=0;
    % Detect Number of blobs/object in the image
    %
    % Choose brightest pixels and look at them
    %
    PickBrightestNPixels=80;
    [sortedValues,sortIndex] = sort(temp3(:),'descend');
    maxIndex = sortIndex(1:PickBrightestNPixels);
    temp5=zeros(512,512);
    temp5(maxIndex)=sortedValues(1:PickBrightestNPixels);
    %
    % Compute blobs
    %
    clear cc
    clear ObjectsSize
    temp6=temp5;
    cc = bwconncomp(temp6);
    %
    % let's pick only the n biggest blobs in the image
    % 
    NumberObjects=5;
    ObjectsSize=zeros(1,length(cc.PixelIdxList));
    TotalDetectedObjects=length(cc.PixelIdxList);
    for iii=1:length(cc.PixelIdxList)
        ObjectsSize(iii)=length(cc.PixelIdxList{iii});
    end
    [sortedValues,sortIndex] = sort(ObjectsSize(:),'descend');
    ObjectsToBeConsidered=min(cc.NumObjects,NumberObjects);
    cc.NumObjects=ObjectsToBeConsidered;
    cc.PixelIdxList=cc.PixelIdxList(sortIndex(1:ObjectsToBeConsidered));
    %
    % Compute centroid of connected components
    %
    S = regionprops(cc,'all');
    %
    % Decide which object to follow
    %
    IsObjectDetected=zeros(1,length(S)); %Object is considered not detected
    AreaObject=zeros(1,length(S));
    IntensityObject=zeros(1,length(S));
    % Let's use our logit function to decide which blob to follow
    distance=zeros(1,length(S));
    % Computing Distance
    for iii=1:length(S)
        y0=S(iii).Centroid(2);
        x0=S(iii).Centroid(1);
        if ~isinf(m)
            distance(iii)=abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/(sqrt(y2-y1)^2+(x2-x1)^2);
        else
            distance(iii)=abs(x0-x1);
        end
        AreaObject(iii)=length(cc.PixelIdxList{iii});
        IntensityObject(iii)=sum(temp5(cc.PixelIdxList{iii}));
    end
    temp7=zeros(512,512);
    % Fitting of the logit function, NB our logit has x1=area,
    % x2=itensity, x3=distance,x4=mean intensity
     X=[AreaObject',IntensityObject',distance',(IntensityObject./AreaObject)'];
    [yhat,dylo,dyhi] = glmval(b,X,'logit',stats);
    for iii=1:length(S)
        if yhat(iii)>=0.4 && yhat(iii)==max(yhat)
            IsObjectDetected(iii)=1;
        end
    end
        %                      X                 Y                Location of the event
        detected=0;
    for iii=1:length(S)
        if IsObjectDetected(iii)==1
        Location=[Location;S(iii).Centroid(1),S(iii).Centroid(2),S(iii).MajorAxisLength,S(iii).Orientation];
        detected=1;
        end
    end
    if detected~=1
        Location=[Location;NaN,NaN,NaN,NaN];
    end
%     
end
output_args=Location;
end