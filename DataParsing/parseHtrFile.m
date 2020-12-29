function Data = parseHtrFile(htrFile)
% Data = parseHtrFile(htrFile)
%
% This routine parses a .htr file into a single data structure for
% subsequent processing
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% htrFile: an absolute or relative path to a htr file
%
%%%%%%%%%%%%%%%%
%%% OUTPUTS %%%%
%%%%%%%%%%%%%%%%
% Data: a structure with fields:
%   nodes: a structure, one per node, with fields
%       Name: the name of the node
%       Parent: the name of the parent node for this node
%       Trans: a [NumberOfFrames x XYZ] matrix of per frame node locations
%       Rot: a [NumberOfFrames x XYZ] matrix of per frame node Euler Angles
%       SF: a [NumberOfFrames x 1] vector of scaling factors
%   params: a structure, with fields corresponding to the header
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com

% check file extension
[fullPath,fileName,ext] = fileparts(htrFile);
if ~strcmp(ext,'.htr')
    error('This is not a .htr file.')
end
if isempty(fullPath)
    fullPath = pwd;
end
params.fullPath = fullPath;
params.fileName = fileName;

% Attempt to open file
fid = fopen(htrFile,'r');
if fid < 0
    error('Invalid .htr file. Check filename.')
end

lineCounter = 0;
while 1
    
    %Scan to file parameters
    tline = fgetl(fid);
    
    %Terminate at the end of file
    if ~ischar(tline)
        fclose(fid);
        error('End of file hit when scanning to parameters. Check file.');
    end
    
    lineCounter = lineCounter + 1;
    
    % Scan to file parameters
    if strfind(tline,'[Header]')
        
        % Advance line past section delimiter
        tline       = fgetl(fid);
        %Terminate at the end of file
        if ~ischar(tline)
            fclose(fid);
            error('End of file hit when parsing parameters. Check file.');
        end
    
        lineCounter = lineCounter + 1;    
        
        %Parse File Parameters
        while isempty(strfind(tline,'[SegmentNames&Hierarchy]'))
            
            % check for comment and split line if necessary
            if strfind(tline,'#')
                splitStr = regexp(tline,'\#','split');
                tline = splitStr{1};
            end
            
            %parameters entries are formatted as {String, String} or as
            %{String, Numeric}
            StrChk = regexp(tline,'[a-z_A-Z]+','match');
            if length(StrChk) == 2
                % {String, String}
                params.(StrChk{1}) = StrChk{2};
            elseif length(StrChk) < 2 && ~isempty(StrChk)
                %{String, Numeric}
                StrChk = regexp(tline,'[a-z_A-Z]+|[0-9.]+','match');
                params.(StrChk{1}) = str2double(StrChk{2});
            else
                fclose(fid);
                error('Error in parsing parameters. More than 2 entries.')
            end
            
            % Advance line
            tline       = fgetl(fid);
            %Terminate at the end of file
            if ~ischar(tline)
                fclose(fid);
                error('End of file hit when parsing parameters. Check file.');
            end
            lineCounter = lineCounter + 1;    
        
        end
    end
    
    % Scan to node connectivity then break to initialize nodes
    if strfind(tline,'[SegmentNames&Hierarchy]')
        break
    end
end

% Initialize nodes
if isfield(params,'NumSegments')
    NumSegments = params.NumSegments;
else
    fclose(fid);
    error('NumSegments is not a parameter in the header file.')
end

nodes = struct('Name',cell(NumSegments,1),'Parent',cell(NumSegments,1),...
               'Trans',cell(NumSegments,1),'Rot',cell(NumSegments,1),...
               'SF',cell(NumSegments,1));

% Make sure end up with the correct number of nodes
nodeCounter = 0;

while 1 
    
    % Scan to node connectivity
    tline = fgetl(fid);
    
    %Terminate at the end of file
    if ~ischar(tline)
        fclose(fid);
        error('End of file hit when scanning to node connectivity. Check file.');
    end
    
    lineCounter = lineCounter + 1;
    
    % Scan to node connectivity
    if strfind(tline,'#CHILD')
        
        % Advance line past section delimiter
        tline       = fgetl(fid);
        %Terminate at the end of file
        if ~ischar(tline)
            fclose(fid);
            error('End of file hit when parsing node connectivity. Check file.');
        end
        
        lineCounter = lineCounter + 1;    
        
        %Parse Node Connectivity
        while isempty(strfind(tline,'[BasePosition]'))
            
            % Make sure end up with the correct number of nodes
            if nodeCounter > NumSegments
                fclose(fid);
                error('Error in parsing node connectivity: too many lines.')
            end
            
            %Log node connectivity
            nodeCounter = nodeCounter + 1;
            StrChk      = regexp(tline,'\w+','match');
            
            if length(StrChk) ~= 2
                fclose(fid);
                error('Error in parsing node connectivity: too many entries.')
            end
            
            nodes(nodeCounter).Name   = StrChk{1};
            nodes(nodeCounter).Parent = StrChk{2};
            
            % Advance line
            tline       = fgetl(fid);
            %Terminate at the end of file
            if ~ischar(tline)
                fclose(fid);
                error('End of file hit when parsing node connectivity. Check file.');
            end
            
            lineCounter = lineCounter + 1;    

        end
    end
    
    % Scan to node connectivity then break to initialize nodes
    if strfind(tline,'[BasePosition]')
        break
    end
end

% Initialize Data Block
if isfield(params,'NumFrames')
    NumFrames = params.NumFrames;
else
    fclose(fid);
    error('NumFrames is not a parameter in the header file.')
end

% The data for each node is formatted as:
% {Tx, Ty, Tz, Rx, Ry, Rz, SF}
DataBlock   = nan(NumFrames,7);
nodeCounter = 0;    % reset node counter when parsing data

while 1

    %Terminate at the end of file
    if ~ischar(tline)
        if nodeCounter < NumSegments
            fclose(fid);
            error('End of file hit when parsing node data. Check file.');
        else
            break
        end
    end
    
    % Scanning to node data
    tline = fgetl(fid);
    
    %Terminate at the end of file
    if ~ischar(tline)
        fclose(fid);
        error('End of file hit when parsing node data. Check file.');
    end
    
    lineCounter = lineCounter + 1;
    
    % Scan to node data section
    if strfind(tline,'#Beginning of Data')
        
        % Advance line past section delimiter
        tline       = fgetl(fid);
        %Terminate at the end of file
        if ~ischar(tline)
            fclose(fid);
            error('End of file hit when parsing node data. Check file.');
        end

        lineCounter = lineCounter + 1;        
            
        %Parse Node Data
        while 1
            
            %Terminate at the end of file
            if ~ischar(tline)
                if nodeCounter < NumSegments
                    fclose(fid);
                    error('End of file hit when parsing node data. Check file.');
                else
                    break
                end
            end
            
            % Parse current Node Name and Increment Node Counter
            StrChk = regexp(tline,'\w+','match');
            % Advance line past sub-section start
            tline       = fgetl(fid);
            %Terminate at the end of file
            if ~ischar(tline)
                fclose(fid);
                error('End of file hit when parsing node data. Check file.');
            end
            % Advance line past sub-section header
            tline       = fgetl(fid);
            %Terminate at the end of file
            if ~ischar(tline)
                fclose(fid);
                error('End of file hit when parsing node data. Check file.');
            end
            lineCounter = lineCounter + 1;
                
            nodeCounter = nodeCounter + 1;
            
            fprintf(1,'Writing out node %d, name: %s.\n',nodeCounter, StrChk{1});
            
            % Advance until the next node, or the file terminates
            while isempty(strfind(tline,'['))
                
                %Parse data for this node and frame
                DataChk = regexp(tline,'[0-9.-]+','match');
                
                % Check for correctness in data
                if length(DataChk) ~= 8
                    fclose(fid);
                    error('Error in parsing node data: formatting error.')
                end
                
                %Log data for this node and frame
                Frame              = str2double(DataChk{1});
                DataBlock(Frame,1) = str2double(DataChk{2});
                DataBlock(Frame,2) = str2double(DataChk{3});
                DataBlock(Frame,3) = str2double(DataChk{4});
                DataBlock(Frame,4) = str2double(DataChk{5});
                DataBlock(Frame,5) = str2double(DataChk{6});
                DataBlock(Frame,6) = str2double(DataChk{7});
                DataBlock(Frame,7) = str2double(DataChk{8});
                
                % Advance line
                tline       = fgetl(fid);
                %Terminate at the end of file
                if ~ischar(tline)
                    if nodeCounter < NumSegments
                        fclose(fid);
                        error('End of file hit when parsing node data. Check file.');
                    else
                        break
                    end
                end
                
                lineCounter = lineCounter + 1;    
                
            end
            
            % Log this Data Block with the appropriate Node
            for ii = 1:NumSegments
                if strcmp(nodes(ii).Name, StrChk{1})
                    nodes(ii).Trans = DataBlock(:,1:3);
                    nodes(ii).Rot   = DataBlock(:,4:6);
                    nodes(ii).SF    = DataBlock(:,7);
                    
                    % Reset DataBlock
                    DataBlock   = nan(NumFrames,7);
                    break
                end
            end
            
        end
    end
    
end

% Compose output
Data.params = params;
Data.nodes  = nodes;

fclose(fid);