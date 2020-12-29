function Map = MapKinect2Simm(mapFile, simData)
% Map = MapKinect2Simm(mapFile, simData)
%
% This routine parses a map file between Kinect2 and SIMM data. It uses
% this information, in conjunction with the simData structure to output a
% N x 2, [KinectJoint SimmJoint], array of correspondences.
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% mapFile: an absolute or relative path to a map file
% simData: a SIMM data structure formatted according to parseTrcFile.m
%
%%%%%%%%%%%%%%%%
%%% OUTPUTS %%%%
%%%%%%%%%%%%%%%%
% Map: a N x 2 array of correspondences between Kinect2 and Simm nodes. The
%   left column is integers indicating the Kinect2 nodes, according to the
%   Kinect2 API, and the right integer index into the simData
%   structure.
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com

% check file extension
[~,~,ext] = fileparts(mapFile);
if ~strcmp(ext,'.txt')
    error('This is not a .txt file and maybe not a map file.')
end

% Attempt to open file
fid = fopen(mapFile,'r');
if fid < 0
    error('Invalid mapfile. Check filename.')
end

if nargin < 2
    error('This routine requires 2 inputs: a map file and a simData structure.');
end

% Kinect2 joints are numbered 1-25
Map         = zeros(25,2);
mapCounter  = 1;

lineCounter = 0;
while 1
    
    %Scan to file parameters
    Msg   = 'End of file hit when scanning to map. Check file.';
    tline = AdvanceLines(fid, 1, Msg);
    
    % Scan to file parameters
    if strfind(tline,'[Kinect2 : SIMM]')
        
        % Advance line past section delimiter
        Msg   = 'End of file hit when parsing map. Check file.';
        tline = AdvanceLines(fid, 1, Msg);
        
        % Kinect2 joints number 1-25
        KinectJointCounter = 1;
        
        % Find Mapping
        while 1
            
            % Mapped joints are delimited by :
            if strfind(tline,':')
                
                splitStr = regexp(tline,'\:','split');
                
                % node label to search for in simData
                chk = regexp(splitStr{2},'[a-z_A-Z_0-9]+','match');
                Map(mapCounter,:) = [KinectJointCounter, ...
                                     find(ismember(simData.labels,chk))];
                
                mapCounter        = mapCounter + 1;
                
            end
            
            % Advance line and check for end of file
            tline       = fgetl(fid);
            %Terminate at the end of file
            if ~ischar(tline)
                fclose(fid);
                Map = Map(1:mapCounter-1,:);
                return
            end
            KinectJointCounter = KinectJointCounter + 1;
        end
    end
end

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