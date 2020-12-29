function SyncInfo = parseAncFile(ancFile)
% Data = parseAncFile(htrFile)
%
% This routine parses a .anc file into a single data structure for
% synchronizing between the Cortex and Kinect2 systems. The .anc file is a
% data log file. The input to this file consists of square pulses sampled
% at the frame rate of the Cortex acquisition system.
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% ancFile: an absolute or relative path to a anc file
%
%%%%%%%%%%%%%%%%
%%% OUTPUTS %%%%
%%%%%%%%%%%%%%%%
% SyncInfo: a structure with fields:
%   syncFrames: integer indices of Cortex to Kinect2 frame correspondences,
%       e.g. 5,10,12,13 means these Cortex frames map to Kinect2 frames
%       1,2,3,4.
%   params: a structure, with fields corresponding to the header
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com

% check file extension
[fullPath,fileName,ext] = fileparts(ancFile);
if ~strcmp(ext,'.anc')
    error('This is not a .htr file.')
end
if isempty(fullPath)
    fullPath = pwd;
end
params.fullPath = fullPath;
params.fileName = fileName;

% Attempt to open file
fid = fopen(ancFile,'r');
if fid < 0
    error('Invalid .htr file. Check filename.')
end

lineCounter = 0;
while 1
    
    %Scan to file parameters
    Msg   = 'End of file hit when scanning to parameters. Check file.';
    tline = AdvanceLines(fid, 1, Msg);
    
    % Scan to file parameters
    if strfind(tline,'File_Type')
        
        % Advance line past section delimiter
        Msg   = 'End of file hit when parsing parameters. Check file.';
        tline = AdvanceLines(fid, 1, Msg);
        
        %Parse File Parameters
        while isempty(strfind(tline,'Range'))
            
            % To find parameters split lines along tab delimiters and find 
            % by knowing that parameters are in key/value pairs.
            
            % The first thing we want to do is to check the file duration
            % so we can pre-allocate memory
            if strfind(tline,'Duration(Sec.)')
                
                splitstr = regexp(tline,'\t','split');
                
                for ii = 1:length(splitstr)
                    if strfind(splitstr{ii},'Duration(Sec.)')
                        params.DurationCortex = str2double(splitstr{ii+1});
                    end
                end
                
            end
            
            % Now let's get our BitDepth and PreciseRate
            if strfind(tline,'BitDepth')
                
                splitstr = regexp(tline,'\t','split');
                
                for ii = 1:length(splitstr)
                    if strfind(splitstr{ii},'BitDepth')
                        params.BitDepth = str2double(splitstr{ii+1});
                    end
                    
                    if strfind(splitstr{ii},'PreciseRate')
                        params.FpsCortex = str2double(splitstr{ii+1});
                    end
                end
                
            end
            
            % Advance line
            Msg   = 'End of file hit when parsing parameters. Check file.';
            tline = AdvanceLines(fid, 1, Msg);
            
        end
    end
    
    % Scan to sync log then break
    if strfind(tline,'Range')
        break
    end
end

% Initialize memory for sync values as the Number of Cortex Frames
syncFrames   = zeros(ceil(params.DurationCortex*params.FpsCortex),1,'int32');
FrameCounter = 0;   % count Cortex frames using the file, discounting the header
syncCounter  = 0;   % counter for our syncFrames buffer

while 1 
    
    % Scan to sync log (This should be the line after 'Range')
    tline = fgetl(fid);
    
    %Terminate at the end of file
    if ~ischar(tline)
        break
    end
    
    lineCounter  = lineCounter + 1;
    FrameCounter = FrameCounter + 1;
    
    % we are looking for positive going values
    if isempty(regexp(tline,'-','once'))
        
        % increment buffer counter and log sync frame
        syncCounter             = syncCounter + 1;
        syncFrames(syncCounter) = FrameCounter;
        
    end
end
% truncate syncFrame buffer to number of detected frames
syncFrames = syncFrames(1:syncCounter);

% Compose output
SyncInfo.params      = params;
SyncInfo.syncFrames  = syncFrames;

fclose(fid);

% NESTED FUNCTIONS
function NextLine = AdvanceLines(fid,Nlines,MsgStr)
    
    advance  = 0;
    while advance < Nlines
        
        NextLine = fgetl(fid);
        
        %Terminate at the end of file
        if ~ischar(NextLine)
            fclose(fid);
            error(MsgStr);
        end
        advance     = advance + 1;
        lineCounter = lineCounter + 1;
    end
end

end