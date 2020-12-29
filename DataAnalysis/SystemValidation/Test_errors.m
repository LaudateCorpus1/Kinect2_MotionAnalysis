% TEST alignment of Kinect2 and SIMM skeletal tracking data
clear
clc

% Test data folder
folder  = 'D:\Kinect2\Data\Data10.10.2014\Kathleen\';
% Test files
trcData     = 'Kathleen_sitting1.mat';
kinect2Data = 'Sitting_Dual.bin';
ancFile     = 'Kathleen_sitting1.anc';

ErrorTitle  = 'Subject3 Sitting';

% map file Kinect2 and SIMM
mapFile = 'D:\Kinect2\DataParsing\MappingKinect2ToSimmJoints.txt';

% Load pre-computed htr file data structure
loadedData  = load(sprintf('%s%s',folder, trcData));
simData     = loadedData.simData;
Nsync       = simData.params.NumFrames;

% Load Kinect2 data
fid         = fopen(sprintf('%s%s',folder,kinect2Data),'r');
fileInfo    = dir(sprintf('%s%s',folder,kinect2Data));
NJoints     = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(NJoints*3*SizeOfFloat);

% Compose Data as [3*NJoints NFrames] matrix and transpose
K2Data      = fread(fid,[3*NJoints NFrames],'float');
K2Data      = K2Data';
fclose(fid);

% Load synchronization file and align frames
SyncInfo    = parseAncFile(sprintf('%s%s',folder,ancFile));

% Parse Map file between kinect2 and SIMM
Map = MapKinect2Simm(mapFile, simData);

% Access aligned frames for Kinect2 and SIMM data
Frames = SyncInfo.syncFrames(SyncInfo.syncFrames < Nsync + 1);
%Kinect2
inds   = [3*(Map(:,1)-1)+1 3*(Map(:,1)-1)+2 3*Map(:,1)]';
inds   = inds(:)';
K2Data = K2Data(1:numel(Frames),inds);
%SIMM
inds     = [3*(Map(:,2)-1)+1 3*(Map(:,2)-1)+2 3*Map(:,2)]';
inds     = inds(:)';
SIMMData = simData.data(Frames,inds);

NFrames = numel(Frames);
errors  = zeros(NFrames, NJoints);
for ii = 1:numel(Frames)
    if ~mod(ii,200)
        fprintf(1,'Processing Frame %d out of %d\n',ii, NFrames);
    end
    
    K2      = reshape(K2Data(ii,:),[3 25]);
    SIMM    = reshape(SIMMData(ii,:),[3 25]);
    
    % estimate least squares solution
    [s,R,T] = estimateSimTrans(K2,SIMM);
    
    % transform Kinect2 data
    K2      = s*R*K2 + repmat(T,[1 25]);
    
    % compute errors
    errors(ii,:) = sqrt(sum((SIMM-K2).^2));
end

labels = cell(size(simData.labels));
for ii = 1:size(simData.labels)
    labels{ii} = simData.labels{Map(ii,2)}(3:end);
end

figure, boxplot(errors/10,'labels',labels,'labelorientation','inline')
    ylabel('Error Distributions (cm)')
    title(sprintf('SIMM versus Kinect2 node location error distributions. File: %s.',ErrorTitle),...
        'FontSize',12)
    set(gcf,'PaperPositionMode','auto')