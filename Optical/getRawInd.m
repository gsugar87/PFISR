function [firstRawIndex,lastRawIndex,ReqRawInd] = getRawInd(BigFN,BytesPerImage,nHeadBytes,ReqFrameInd)
% Michael Hirsch 2014
% gets first and last raw indicies from a big .DMCdata file

Nmetadata = nHeadBytes/2; %number of 16-bit words
fid = fopen(BigFN,'r');
%% get first raw frame index
fseek(fid,BytesPerImage,'bof');
firstRawIndex = meta2rawInd(fid,Nmetadata);
%% get last raw frame index
fseek(fid, -nHeadBytes, 'eof');
lastRawIndex = meta2rawInd(fid,Nmetadata);
%% take requests
if nargin > 3 && nargout > 2
    nFrameReq = length(ReqFrameInd);
    
    ReqRawInd = zeros(nFrameReq,1,'double'); %preallocate, double float for numerical ops compatibility
    jFrm = 0;
    for iFrm = ReqFrameInd
        jFrm = jFrm + 1;
        
        currByte = (iFrm - 0) * BytesPerImage;  %remember, we want the header, not the image
        fseek(fid,currByte,'bof');
        ReqRawInd(jFrm) = meta2rawInd(fid,Nmetadata);
    end %for
end %if
%%
fclose(fid);
end %function

function rawInd = meta2rawInd(fid,Nmetadata)
    metadata = fread(fid,Nmetadata,'uint16=>uint16',0,'l');
    
    %typecast those 16-bit metadata words into frame #'s
    rawInd = typecast([metadata(2) metadata(1)],'uint32');
    rawInd = double(rawInd); %to do expected math operations
end %function