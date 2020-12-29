function NormalizeKinect2Data(kinect2BinFile)
%NormalizeKinect2Data : 
% This file takes in a binary file corresponding
% to data acquired via the kinect2Acq system and writes out a mp4 video file
% for inspecting the quality of the data
%
%%% INPUTS %%%
%
% kinect2BinFile : a string path to the kinect2 binary file
% 
%%% OUTPUTS %%%
%
% a .mat file in the same location, named according to the kinect2
% binary file, containing the normlaized kinematic model data
%
% author: John D. Long II, PhD, contact: jlong29@gmail.com


mfname = mfilename;

% check if file exists
fileInfo = dir(kinect2BinFile);
if isempty(fileInfo)
    error(['%s: The file ''%s'' was not found on this path.'...
            ' Did you suppply the correct path? Is this file missing?'], mfname, kinect2BinFile);
end

% Check file extension
[~, filename, ext] = fileparts(kinect2BinFile);

if ~strcmp(ext, '.bin')
    error('%s: The file ''%s'' is not a binary file.', mfname, kinect2BinFile);
end

% Load Kinect2 data
fid         = fopen(kinect2BinFile,'r');
Nnodes      = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(Nnodes*3*SizeOfFloat);

% Read in data
K2Data      = fread(fid,[3*Nnodes NFrames],'float');
fclose(fid);

% check correctness of file
if (NFrames - floor(NFrames)) > 0
    error(['%s: Incorrect File: The computed number of frames in ''%s'' is '...
           'not an integer value.'], mfname, kinect2BinFile);
end
fprintf(1,'Number of frames in Kinect2 file %s: %d.\n',kinect2BinFile, NFrames);

% Create kinematic chain dependency graph according to documentation:
% "Kinect for Windows SDK 2.0" It is zero based, so I'll correct.
KmSegments    = cell(7,1); % there are 7 segments to be accounted for
KmSegments{1} = [0 12 13 14 15]' + 1;   % base of spine down left leg
KmSegments{2} = [0 16 17 18 19]' + 1;   % base of spine down right leg
KmSegments{3} = [0 1 20 2 3]' + 1;      % base of spine to tip of head
KmSegments{4} = [20 4 5 6 7 22]' + 1;   % left shoulder through thumb
KmSegments{5} = [7 21]' + 1;            % left hand to finger tip
KmSegments{6} = [20 8 9 10 11 24]' + 1; % right shoulder through thumb
KmSegments{7} = [11 23]' + 1;           % right hand to finger tip

% There is 3 special cases to be handled. When making adjustments from the
% base of the spine through the spine shoulder, the left and right arms
% must inherit these translations. Likewise, when scaling from the
% shoulders through the thumbs, the finger tips must be adjusted
% accordingly.
Case1 = [KmSegments{4}(2:end);
         KmSegments{5}(2:end);
         KmSegments{6}(2:end);
         KmSegments{7}(2:end)'];

Case2 = 22;
Case3 = 24;

% Following the segments detailed above, specify all nodes pairs in an edge
Connections = nan(100,2);
ElemCounter = 1;    
for ii = 1:length(KmSegments)
    tmp         = KmSegments{ii};
    Connections(ElemCounter:ElemCounter + length(tmp)-2,:) = [tmp(1:end-1) tmp(2:end)];
    ElemCounter = ElemCounter + length(tmp) - 1;
end
ElemCounter = ElemCounter - 1;
Connections = Connections(1:ElemCounter,:);

EdgeCounter = 1;
K2Scaled    = nan(size(K2Data));
% Access each frame
for ii = 1:NFrames
    
    % Access this frame in R3
    K2             = reshape(K2Data(:,ii),[3 Nnodes]);
    
    % Translate base of spine to origin
    TranslateKinect2Data();
    
    % Rotate kinect2 data such that the spine runs along the y-axis and the
    % line from the base of the spine through the right hip runs along the
    % x-axis
    RotateKinect2Data();
    
    % Scale Kinect2 data such that all segments have unit length
    ScaleKinect2Data();
    
    % Log Scaled Data
    K2Scaled(:,ii) = K2(:);
end

fprintf(1,'%s: Writing out binary file ''%s''...',mfname, [filename '_normalized' ext]);
% Write out normalized data to file
fid = fopen([filename '_normalized' ext],'wb');
% Compose Data as [3*Nnodes NFrames] matrix and transpose
chk = fwrite(fid,K2Scaled,'float');
fprintf(1,'done.\n');
if chk ~= (3*Nnodes*NFrames)
    fclose(fid);
    error(['%s: File Write Error: The file %s was not normalized because '...
           'only %d out of %d frames were written to file.'], mfname, kinect2BinFile, chk, NFrames);
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%
%%% NESTED Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%
function TranslateKinect2Data()
    K2 = K2 - repmat(K2(:,1),[1 Nnodes]);
end

function RotateKinect2Data()
    % 1) Align the axis from the base of the spine to the right hip with (1,0,0)
    pose   = K2;     % debugging copy
    axis1  = [1 0 0]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Rotation to axis1 %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    ModRef   = pose(:,17) - pose(:,1);
    ModRef   = ModRef./norm(ModRef);
    crossCM  = cross(ModRef,axis1);
    n        = crossCM./norm(crossCM);  % plane normal (rotation axis)
    theta    = real(asin(n'*crossCM));
    T        = twist_matrix(n,[0 0 0]',theta);
    
    % Transform pose to align with axis1
    pose = T*cat(1,pose,ones(1,Nnodes));
    pose = pose(1:3,:);
    
    % 3) Align axis from base of spine to lower back with (0,1,0)
    
    axis2  = [0 1 0]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Rotation to axis2 %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Looking down the axis from the right hip to the base of the spine, align
    % the segment from the base of the spine to the lower back with the y-axis
    
    % using MIRM (Murray and Sastry p. 100)
    w      = axis1;
    u_hat  = pose(:,2) - pose(:,1); u_hat(1) = 0;
    v_hat  = axis2;
    
    theta  = atan2(w'*cross(u_hat,v_hat),u_hat'*v_hat);
    % Compose Transform
    U      = twist_matrix(w,zeros(3,1),theta);
    % Transform Points
    pose = U*cat(1,pose,ones(1,Nnodes));
    
    % Lastly, a bit of rotation in the xz-plane to make sure the spine
    % segment aligns with the y-axis. This will misalign the right hip a
    % bit from the x-axis.
    ModRef = pose(1:3,2) - pose(1:3,1);
    ModRef = ModRef./norm(ModRef);
    theta  = acos([0 1 0]*ModRef);
    V      = twist_matrix([0 0 1]',zeros(3,1), -theta);
    % Transform Points
    pose = V*pose;
    
    % copy back final result
    K2   = pose(1:3,:);
end

function ScaleKinect2Data()
    SegmentLengths = ...
    sqrt(sum((K2(:,Connections(:,2)) - K2(:,Connections(:,1))).^2));

    % Access all kinematic chain segments and make adjustments
    for jj = 1:length(KmSegments)

        adjust = KmSegments{jj};

        for kk = 2:length(adjust)
            % Access current direction vector between 2 nodes
            D                    = K2(:,adjust(kk))-K2(:,adjust(kk - 1));
            % Adjust this edge to be unit length and translate all downstream nodes
            offset               = (1/norm(SegmentLengths(EdgeCounter))).*D - D;
            K2(:,adjust(kk:end)) = K2(:,adjust(kk:end)) + repmat(offset,[1 length(adjust(kk:end))]);
            % Increment EdgeCounter
            EdgeCounter          = EdgeCounter + 1;

            % handle forks

            % Case1 : Base of spine through shoulder spine
            if jj == 3
                if adjust(kk) == 2 || adjust(kk) == 21
                    K2(:,Case1) = K2(:,Case1) + repmat(offset,[1 length(Case1)]);
                end 
            end
            % Case2 : left shoulder through thumb
            if jj == 4
                if adjust(kk) ~= 23
                    K2(:,Case2) = K2(:,Case2) + repmat(offset,[1 length(Case2)]);
                end 
            end
            % Case3 : right shoulder through thumb
            if jj == 6
                if adjust(kk) ~= 25
                    K2(:,Case3) = K2(:,Case3) + repmat(offset,[1 length(Case3)]);
                end 
            end
        end
    end
    %reset EdgeCounter
    EdgeCounter = 1;
end

end %EOF