function [data] = SEV2mat(SEV_DIR, varargin)
%SEV2MAT  TDT SEV file format extraction.
%   data = SEV2mat(SEV_DIR), where SEV_DIR is a string, retrieves
%   all sev data from specified directory in struct format. SEV files
%   are generated by an RS4 Data Streamer, or by setting the Unique
%   Channel Files option in Stream_Store_MC or Stream_Store_MC2 macro
%   to Yes.
%
%   data    contains all continuous data (sampling rate and raw data)
%
%   data = SEV2mat(SEV_DIR,'parameter',value,...)
%
%   'parameter', value pairs
%      'CHANNEL'    integer, returns the sev data from specified channel
%                       only (default = 0 for all channels)
%      'JUSTNAMES'  boolean, retrieve only the valid event names
%      'EVENTNAME'  string, specific event name to retrieve data from
%      'VERBOSE'    boolean, set to false to disable console output
%      'DEVICE'     string, connect to specific RS4 device.  DEVICE can be
%                       the IP address or NetBIOS name of RS4-device 
%                       (e.g. RS4-41001).  Requires TANK and BLOCK
%                       parameters
%      'TANK'       string, tank on RS4 to retrieve data from. Requires
%                       DEVICE and BLOCK parameters
%      'BLOCK'      string, block on RS4 to retrieve data from. Requires
%                       DEVICE and TANK parameters

% defaults
CHANNEL   = 0;
EVENTNAME = '';
DEVICE    = '';
TANK      = '';
BLOCK     = '';
VERBOSE   = 1;
JUSTNAMES = 0;

% parse varargin
for i = 1:2:length(varargin)
    eval([upper(varargin{i}) '=varargin{i+1};']);
end

if any([~isempty(DEVICE) ~isempty(TANK) ~isempty(BLOCK)])
    if any([isempty(DEVICE) isempty(TANK) isempty(BLOCK)])
        error('DEVICE, TANK and BLOCK must all be specified');
    else
        SEV_DIR = sprintf('\\\\%s\\data\\%s\\%s\\', DEVICE, TANK, BLOCK);
    end
end

data = [];

ALLOWED_FORMATS = {'single','int32','int16','int8','double','int64'};

eventNames = {};
eventNameCount = 0;

if strcmp(SEV_DIR(end), filesep) == 0
    SEV_DIR = [SEV_DIR filesep];
end

file_list = dir([SEV_DIR '*.sev']);
if length(file_list) < 1
    warning(['no sev files found in ' SEV_DIR])
    return
end

for i = 1:length(file_list)
    path = [SEV_DIR file_list(i).name];
    
    % open file
    fid = fopen(path, 'rb');
    
    if fid < 0
        warning([path ' not opened'])
        return
    end
    
    % create and fill streamHeader struct
    streamHeader = [];
    
    streamHeader.fileSizeBytes   = fread(fid,1,'uint64');
    streamHeader.fileType        = char(fread(fid,3,'char')');
    streamHeader.fileVersion     = fread(fid,1,'char');
    
    if streamHeader.fileVersion < 3
        
        % event name of stream
        if streamHeader.fileVersion == 2 
            streamHeader.eventName  = char(fread(fid,4,'char')');
        else
            streamHeader.eventName  = fliplr(char(fread(fid,4,'char')'));
        end
        
        % current channel of stream
        streamHeader.channelNum        = fread(fid, 1, 'uint16');
        % total number of channels in the stream
        streamHeader.totalNumChannels  = fread(fid, 1, 'uint16');
        % number of bytes per sample
        streamHeader.sampleWidthBytes  = fread(fid, 1, 'uint16');
        reserved                 = fread(fid, 1, 'uint16');
        
        % data format of stream in lower four bits
        streamHeader.dForm      = ALLOWED_FORMATS{bitand(fread(fid, 1, 'uint8'),7)+1};
        
        % used to compute actual sampling rate
        streamHeader.decimate   = fread(fid, 1, 'uint8');
        streamHeader.rate       = fread(fid, 1, 'uint16');
        
        % reserved tags
        reserved = fread(fid, 1, 'uint64');
        reserved = fread(fid, 2, 'uint16');
       
    else
        error(['unknown version ' num2str(streamHeader.fileVersion)]);
    end
    
    if streamHeader.fileVersion > 0
        % determine data sampling rate
        streamHeader.fs = 2^(streamHeader.rate)*25000000/2^12/streamHeader.decimate;
        % handle multiple data streams in one folder
        exists = isfield(data, streamHeader.eventName);
    else
        streamHeader.dForm = 'single';
        streamHeader.fs = 0;
        s = regexp(file_list(i).name, '_', 'split');
        streamHeader.eventName = s{end-1};
        streamHeader.channelNum = str2double(regexp(s{end},  '\d+', 'match'));
        warning('%s has empty header; assuming %s ch %d format %s and fs = %.2f\nupgrade to OpenEx v2.18 or above\n', ...
            file_list(i).name, streamHeader.eventName, ...
            streamHeader.channelNum, streamHeader.dForm, 24414.0625);
        
        exists = 1;
        %data.(streamHeader.eventName).fs = streamHeader.fs;
        data.(streamHeader.eventName).fs = 24414.0625;
    end

    % skip if this isn't exactly what we're looking for
    if ~strcmp(EVENTNAME, '') && ~strcmp(EVENTNAME, streamHeader.eventName)
        fclose(fid);
        continue
    end
    
    % read rest of file into data array as correct format
    
    if JUSTNAMES
        bFoundIt = 0;
        for name = 1:length(eventNames)
            bFoundIt = strcmp(eventNames{name}, streamHeader.eventName);
        end
        if bFoundIt == 0
            eventNameCount = eventNameCount + 1;
            eventNames{eventNameCount} = streamHeader.eventName;
        end
    else
        data.(streamHeader.eventName).name = streamHeader.eventName;
        if CHANNEL > 0
            if streamHeader.channelNum == CHANNEL
                temp_data = fread(fid, inf, ['*' streamHeader.dForm])';
                data.(streamHeader.eventName).data = temp_data;
                data.(streamHeader.eventName).fs = streamHeader.fs;
            end
        else
            if exists ~= 1
                %preallocate data array
                temp_data = fread(fid, inf, ['*' streamHeader.dForm])';
                total_samples = length(temp_data);
                func = str2func(streamHeader.dForm);
                data.(streamHeader.eventName).data = func(zeros(streamHeader.totalNumChannels,total_samples));
                data.(streamHeader.eventName).data(streamHeader.channelNum,:) = temp_data;
                data.(streamHeader.eventName).fs = streamHeader.fs;
            else
                data.(streamHeader.eventName).data(streamHeader.channelNum,:) = fread(fid, inf, ['*' streamHeader.dForm])';
            end
        end
    end
    
    % close file
    fclose(fid);
    %path
    if VERBOSE
        streamHeader
    end
    
    %if streamHeader.fileVersion > 0
    %    % verify streamHeader is 40 bytes
    %    dataSize = length(streamData) * streamHeader.sampleWidthBytes;
    %    streamHeaderSizeBytes = streamHeader.fileSizeBytes - dataSize;
    %    if streamHeaderSizeBytes ~= 40
    %        warning('streamHeader Size Mismatch -- %d bytes vs 40 bytes', streamHeaderSizeBytes);
    %    end
    %end
end

if JUSTNAMES, data = eventNames; end
end

