function [non_zero_elements_meteora,m,q]=FindPixelsOfInterestForMeteor(xpix1,ypix1,xpix2,ypix2)

initialLocation=[xpix1 ypix1 0 xpix2 ypix2];


m=(initialLocation(2)-initialLocation(5))/(initialLocation(1)-initialLocation(4));
q=(initialLocation(2)-initialLocation(5))/(initialLocation(1)-initialLocation(4))*(-initialLocation(4))+initialLocation(5);
x=[1:1:512];
meteora_back_mean_diag=zeros(512,512);
if ~isinf(m) && m~=0
    for qprime=[-(abs(round(0.05*q))+abs(1/m)):1:abs(round(0.05*q))+abs(1/m)]
        x=[1:1:512];
        y=m*x+q+qprime;
        numerazione=find(y>=1&y<=512);
        x=x(numerazione);
        y=round(y(numerazione));
        for lll=1:length(x)
            meteora_back_mean_diag(y(lll),x(lll))=1;
        end
    end
elseif isinf(m)
    for zzz=[xpix1-12:1:xpix1+12]
        if zzz>0 && zzz<513
            meteora_back_mean_diag(:,zzz)=1;
        end
    end
elseif m==0;
    y=[q-12:q+12];
    y=y(y>0);
    y=y(y<513);
    for zzz=1:length(y)
        meteora_back_mean_diag(y(zzz),:)=1;
    end
end

 non_zero_elements_meteora=find(meteora_back_mean_diag);
%     figure
%     imagesc(meteora_back_mean_diag);
%     pause


if length(non_zero_elements_meteora)<=5000
%
%
%
%   OLD METHOD
% 

    meteora_back_mean_diag=zeros(512,512);
    for eee=-25:1:25
        for qqq=1:length(x)
            www=y(qqq)+eee;
                if www>=1 & www<=512
                    meteora_back_mean_diag(www,x(qqq))=1;
                end
        end
    end
    % figure
    % imagesc(meteora_back_mean_diag);
    % pause

    y=[1:1:512];
    x=(y-q)/m;
    numerazione=find(x>=1&x<=512);
    x=round(x(numerazione));
    y=(y(numerazione));
    for eee=-25:1:25
        for qqq=1:length(x)
            www=x(qqq)+eee;
                if www>=1 & www<=512
                    meteora_back_mean_diag(y(qqq),www)=1;
                end
        end
    end
    non_zero_elements_meteora=find(meteora_back_mean_diag);
    % figure
    % imagesc(meteora_back_mean_diag);
    % pause
end

end