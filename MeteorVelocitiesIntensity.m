function [ Intensity] = MeteorVelocitiesIntensity( Location,data_event,data_baseline,var_baseline)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

IntensityForeground=zeros(size(data_event,3),1);
IntensityBackground=IntensityForeground;
tempvar=var_baseline;
% tempstd=var(double(dataUltra),0,3);
dataUltraBackground=data_baseline;
dataUltra=data_event;
% end
% for boxsize=6;
   boxsize=10;
	tempcounter=0;
    for iii=1:size(dataUltra,3)%86
%         iii
%         size(dataUltra)
        if ~isnan(Location(iii,1)) && ~isnan(Location(iii,2)) && (Location(iii,1))>=1 && (Location(iii,1))<=512 ...
                && (Location(iii,2))>=1 && (Location(iii,2))<=512
            tempcounter=tempcounter+1;
            boundsy=round([Location(iii,2)-boxsize,Location(iii,2)+boxsize]);
            boundsx=round([Location(iii,1)-boxsize,Location(iii,1)+boxsize]);
            if boundsx(1)<1
                boundsx(1)=1;
            end
            if boundsx(1)>512;
                boundsx(1)=512;
            end
            if boundsx(end)<1
                boundsx(end)=1;
            end
            if boundsx(end)>512;
                boundsx(end)=512;
            end
            if boundsy(1)<1
                boundsy(1)=1;
            end
            if boundsy(1)>512;
                boundsy(1)=512;
            end
            if boundsy(end)<1
                boundsy(end)=1;
            end
            if boundsy(end)>512;
                boundsy(end)=512;
            end

            if tempcounter<2
                %brightness of meteor
               IntensityForeground(iii)=sum(sum(dataUltra(boundsy(1):boundsy(end),...
               boundsx(1):boundsx(end),iii)));
            else
                %Change in intensity due to meteor, subtracting what is
                %being produced at the previous time step
            IntensityForeground(iii)=sum(sum(dataUltra(boundsy(1):boundsy(end),...
               boundsx(1):boundsx(end),iii)));
%             IntensityForegroundDueToOtherLight(iii)=(sum(sum(dataUltra(boundsy(1):boundsy(end),...
%             boundsx(1):boundsx(end),iii-1)))-sum(sum(dataUltraBackground(boundsy(1):boundsy(end),...
%             boundsx(1):boundsx(end)))));
            end
            % Brightness of background
            IntensityBackground(iii)=sum(sum(dataUltraBackground(boundsy(1):boundsy(end),...
                boundsx(1):boundsx(end))));
            Intensity(iii)=IntensityForeground(iii)-IntensityBackground(iii);
%             MaskBack=dataUltraBackground(boundsy(1):boundsy(end),...
%                 boundsx(1):boundsx(end));
%             Mask=dataUltra(boundsy(1):boundsy(end),boundsx(1):boundsx(end),iii);
%             figure(1)
%             imagesc(double(Mask)-MaskBack)
%             figure(2)
%             BW = edge(double(Mask)-MaskBack)
%             imagesc(BW)
%             pause
                       
        end
    end
% end
end

