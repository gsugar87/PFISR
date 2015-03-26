
%
%
%       Let's see if I can recognize any meteor
%

% Playing first 100 frames of video
close all
clear all
maxzzz=3500;
numFramesAvg = 3500;
RBrEv=zeros(round((maxzzz+numFramesAvg)*.3),4);
rbrcounter=1;
VBrEv=zeros(round((maxzzz+numFramesAvg)*.2),4);
vbrevcounter=1;
EBrEv=zeros(round((maxzzz+numFramesAvg)*.2),4);
ebrevcounter=1;

for zzz=numFramesAvg:numFramesAvg:maxzzz
    [data,~,tUT] = rawDMCreader('D:\2014-04-01\2014-04-01T06-00-CamSer1387.DMCdata',512,512,1,1,zzz:zzz+numFramesAvg,1,[1000, 1100],'auto','auto');
    [data,~,tUT] = rawDMCreader('D:\2014-04-01\2014-04-01T06-00-CamSer1387.DMCdata',512,512,1,1,a(zzz,3),0,[1000, 1100],'auto','auto');

    % % % 1 tiny metor discovered on 01/04/2014 at frame 6036 (between
    % 6030-6040)
    dataprova=single(data);
    provamean=mean(data,3);
    provastd=std(dataprova,0,3);
%     provapercentile=prctile(dataprova,99,3);

    % % bright events
    % BrEv=[];
    % for i=1:512
    %     for j=1:512
    %         % dummy variable where I allocate where in a pixel I have significant 
    %         % difference from standard deviation
    %         dummy=find(dataprova(i,j,:)>(provamean(i,j)+5*provastd(i,j)));
    %         if isempty(dummy)==0
    %             % the idea behind recording the size is related to the
    %             % brightness of an event, if event is really bright --> lux
    %             % spread around
    %             BrEv=[BrEv;i,j,dummy(1)+zzz,size(dummy,1)];
    %         end
    %     end
    % end
    

%    really bright event
    for i=1:512
        for j=1:512
            dummy=find(dataprova(i,j,:)>(provamean(i,j)+7*provastd(i,j)));
            if isempty(dummy)==0
                RBrEv(rbrcounter,:)=[i,j,dummy(1)+zzz-1,size(dummy,1)];
                rbrcounter=rbrcounter+1;
            end
        end
    end
    [~,b]=sort(RBrEv(:,3));
    RBrEv1=RBrEv(b,:);
    % % % very bright event
    for i=1:512
        for j=1:512
            dummy=find(dataprova(i,j,:)>(provamean(i,j)+9*provastd(i,j)));
            if isempty(dummy)==0
                VBrEv(vbrevcounter,:)=[i,j,dummy(1)+zzz-1,size(dummy,1)];
                vbrevcounter=vbrevcounter+1;
            end
        end
    end
    [a,b]=sort(VBrEv(:,3));
    VBrEv1=VBrEv(b,:);
    % % % extremely bright event
    for i=1:512
        for j=1:512
            dummy=find(dataprova(i,j,:)>(provamean(i,j)+15*provastd(i,j)));
            if isempty(dummy)==0
                EBrEv(ebrevcounter,:)=[i,j,dummy(1)+zzz-1,size(dummy,1)];
                ebrevcounter=ebrevcounter+1;
            end
        end
    end
    [a,b]=sort(EBrEv(:,3));
    EBrEv1=EBrEv(b,:);
% 
% %     for iii=1:length(EBrEv)
% %         [data,~,tUT] = rawDMCreaderLorenzoLimonta('/pdata/Alaska2014Data/2014-04-01/2014-04-01T07-01-CamSer7196.DMCdata',512,512,1,1,EBrEv(iii,3)-10:EBrEv(iii,3)+10,2,[50,1000],'auto','auto',EBrEv(iii,1),EBrEv(iii,2));
% %     end
end
save('RBrEV1.mat','RBrEv02042014')
save('VBrEV1.mat','VBrEv')
save('EBrEV1.mat','EBrEv') 
% zzz=17600-3500
% while zzz>3000
%     zzz=zzz+3500;
%  [data,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-04-01/2014-04-01T07-01-CamSer7196.DMCdata',512,512,1,1,zzz:(zzz+3500),0.01,[100,1000],'auto','auto');
% end

jjj=1;
for iii=2:length(RBrEv1)
    if RBrEv1(iii,3)~=RBrEv1(iii-1,3)
        RBrEv2(jjj,:)=RBrEv1(iii,:);
        jjj=jjj+1;
    end
end


jjj=1;
for iii=2:length(VBrEv1)
    if VBrEv1(iii,3)~=VBrEv1(iii-1,3)
        VBrEv2(jjj,:)=VBrEv1(iii,:);
        jjj=jjj+1;
    end
end

jjj=1;
for iii=2:length(EBrEv1)
    if EBrEv1(iii,3)~=EBrEv1(iii-1,3)
        EBrEv2(jjj,:)=EBrEv1(iii,:);
        jjj=jjj+1;
    end
end

clear RBrEv3
zzz=1;
jjj=1;
a=23;
while zzz<=length(RBrEv2)-a
    if RBrEv2(zzz+a,3)>=(RBrEv2(zzz,3)+a+12)
        RBrEv3(jjj,:)=RBrEv2(zzz,:);
        jjj=jjj+1;
        zzz=zzz+1;
    else
        zzz=zzz+a;
    end
end

jjj=1;
for iii=2:length(RBrEv3)
    if RBrEv3(iii,3)==RBrEv3(iii-1,3)+1
        RBrEv4(jjj,:)=RBrEv3(iii,:);
        jjj=jjj+1;
    end
end



% % 
% for zzz=1:length(RBrEv4)
%     zzz
%     [data,~,tUT] = rawDMCreaderLorenzoLimonta('/media/ExtBook/Ultra/2014-04-01/2014-04-01T07-01-CamSer7196.DMCdata',512,512,1,1,RBrEv4(zzz,3):(RBrEv4(zzz,3)+1),0.01,[100,1000],'auto','auto',RBrEv4(zzz,2),RBrEv4(zzz,1));
% end

% 
% 
% %%%% METEORS
% % 20484 TO 20490
% % 56764
% 
% 
% clear all
% for zzz=160100:20:160100+20
%     [data,~,tUT] = rawDMCreader('/media/ExtBook/Ultra/2014-04-01/2014-04-01T07-01-CamSer7196.DMCdata',512,512,1,1,zzz:zzz+20,0.01,[50, 800],'auto','auto');
% 
%     % % % 1 tiny metor discovered on 01/04/2014 at frame 6036 (between
%     % 6030-6040)
%     dataprova=single(data);
% end
% a=sum(dataprova,3);
% figure
% imagesc(a)

zzz=1;
while (zzz>=1) && (zzz<=length(RBrEv4))
    [data,~,tUT] = rawDMCreaderLorenzoLimonta('/media/ExtBook/Ultra/2014-04-01/2014-04-01T07-01-CamSer7196.DMCdata',512,512,1,1,RBrEv4(zzz,3):(RBrEv4(zzz,3)+1),0.01,[100,1000],'auto','auto',RBrEv4(zzz,2),RBrEv4(zzz,1));
    zzz=zzz+1;
end
