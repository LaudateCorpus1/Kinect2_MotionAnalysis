function TrcData = parseTrcFile(trcFile)
% Data = parseTrcFile(trcFile)
%
% This routine parses a .trc file into a single data structure for
% subsequent processing
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% trcFile: an absolute or relative path to a trc file
%
%%%%%%%%%%%%%%%%
%%% OUTPUTS %%%%
%%%%%%%%%%%%%%%%
% TrcData: a structure with fields:
%   labels: a cell array of strings containing the node labels in the order
%       in which they appear in the data field
%   data: a matrix of dimensions [NumberOfFrames 3*NumberOfLabels], with
%       each row containing the xyz coords for each node listed from left
%       to right according to the order in the labels field.
%   params: a structure, with fields corresponding to the header
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com

% check file extension
[fullPath,fileName,ext] = fileparts(trcFile);
if ~strcmp(ext,'.trc')
    error('This is not a .trc file.')
end
if isempty(fullPath)
    fullPath = pwd;
end
params.fullPath = fullPath;
params.fileName = fileName;

% Attempt to open file
fid = fopen(trcFile,'r');
if fid < 0
    error('Invalid .trc file. Check filename.')
end

lineCounter = 0;
while lineCounter < 3
    
    %Scan Header
    tline = AdvanceLines(fid, 1, 'End of file hit when scanning Header. Check file.');
    
    %Parse Header: the first 3 lines are the header. This is very explicit
    %and frail
    
    % Line 1 contains the file name
    if lineCounter == 1
        
        splitStr    = regexp(tline,'\t','split');
        params.File = splitStr{4};
        
    end
    
    %Line 2 contains the keys and line 3 contains the associated values
    if lineCounter == 2
        
        keys = regexp(tline,'\t','split');
         
        % Advance to the next line to get the values
        tline = AdvanceLines(fid, 1, 'End of file hit when parameters key/values. Check file.');
        
        values = regexp(tline,'\t','split');
        
        if length(keys) ~= length(values)
            error('The number of key/values pairs in the header don''t match. Check file.');
        end
        
        for ii = 1:length(keys)
            % check for numeric data
            chk = regexp(values{ii},'[0-9.]+','match');
            if isempty(chk)
                params.(keys{ii}) = values{ii};
            else
                params.(keys{ii}) = str2double(chk);
            end
        end
    end    
end

% Advance to the next line to Access all node labels
tline = AdvanceLines(fid, 1, 'End of file hit when parameters key/values. Check file.');

% Split up the string by tabs, search for our nodes by tag, while noting
% exclusions
tag        = 'V_';
Exclude    = {'V_MidASIS', 'V_MidHip','V_Pelvic_Origin'};
pattern    = sprintf('(%s)[.a-z_A-Z]+',tag);

% The first 3 entries in the .trc file is separated by \t while all other
% labels are separated by \t\t\t. This is annnoying. Moreover, Frame# and
% Time require 1 float, while all the other labels require 3 floats.
NodeLabels = regexp(tline,'\t\t\t','split');

firstbit   = regexp(NodeLabels{1},'\t','split');

NodeLabels = [firstbit(3) NodeLabels(2:end)];

labels      = cell(100,1);
inds        = nan(100,1);
nodeCounter = 0;
for ii = 1:length(NodeLabels)
    % check for tag
    chk = regexp(NodeLabels{ii},pattern,'match');
    if ~isempty(chk)
        % check exclusions
        if ismember(chk, Exclude)
            continue
        else
            nodeCounter = nodeCounter + 1;
            labels{nodeCounter} = chk{1};
            inds(nodeCounter)   = ii;
        end
    end
end
% truncate to correct size
labels = labels(1:nodeCounter);
inds   = inds(1:nodeCounter);

% Node Data is accessed at xyz coords with +2 offset for Frame# and Time
inds   = reshape([3*inds 3*inds+1 3*inds+2]',[1 3*length(inds)]);
%%%%%%%%%%%%%%%%%%
%%% Parse Data %%%
%%%%%%%%%%%%%%%%%%
% skip lines 5 and 6
Msg    = 'End of file hit when advancing to data. Check file.';
AdvanceLines(fid, 2, Msg);

% Initialize Data matrix output
Data   = nan(params.NumFrames, 3*length(labels));

Msg    = 'End of file hit when scanning to data. Check file.';
FrameCounter = 0;
while FrameCounter < params.NumFrames
    
    tline        = AdvanceLines(fid, 1, Msg);
    FrameCounter = FrameCounter + 1;
    
    if ~mod(FrameCounter,200)
        fprintf(1,'Processing Frame %d out of %d\n',FrameCounter, params.NumFrames);
    end
    
    %Convert line string to numeric data
    Frame        = cellfun(@str2double, regexp(tline,'\t','split'));
    
    % Log Frame Data
    Data(FrameCounter,:) = Frame(inds);
end

% Compose output
TrcData.labels = labels;
TrcData.data   = Data;
TrcData.params = params;

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

end %EOF