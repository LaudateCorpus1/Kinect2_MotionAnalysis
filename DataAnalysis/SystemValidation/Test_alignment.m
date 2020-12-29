% TEST alignment of Kinect2 and SIMM skeletal tracking data
clear
clc

% test data folder
folder  = 'D:\Kinect2\Data\Data10.10.2014\Kathleen\';
% test files
trcData     = 'Kathleen_sitting1.mat';
kinect2Data = 'Sitting_Dual.bin';
ancFile     = 'Kathleen_sitting1.anc';

VideoTitle  = 'KathleenSitting10102014';
fps         = 60/2; % num is max fps and div is slowing factor

% map file Kinect2 and SIMM
mapFile = 'D:\Kinect2\DataParsing\MappingKinect2ToSimmJoints.txt';

% Either Load pre-computed trc file data structure, or compute it
if isempty(dir(trcData))
    fprintf(1,'Computing SIMMdata from trc file...\n');
    simData = parseTrcFile([trcData(1:end-4) '.trc']);
    save([folder trcData], 'simData');
else
    fprintf(1,'Loading stored SIMdata structure for this subject.\n');
    loadedData  = load(sprintf('%s%s',folder, trcData));
    simData     = loadedData.simData;
end
Nsync       = simData.params.NumFrames;

% Load Kinect2 data
fid         = fopen(sprintf('%s%s',folder,kinect2Data),'r');
fileInfo    = dir(sprintf('%s%s',folder,kinect2Data));
NJoints     = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(NJoints*3*SizeOfFloat);
fprintf(1,'Number of frames in Kinect2 file: %d.\n',NFrames);

% Compose Data as [3*NJoints NFrames] matrix and transpose
K2Data      = fread(fid,[3*NJoints NFrames],'float');
K2Data      = K2Data';
fclose(fid);

% Load synchronization file and align frames
SyncInfo    = parseAncFile(sprintf('%s%s',folder,ancFile));
fprintf(1,'Number of synched frames in .anc file: %d.\n',length(SyncInfo.syncFrames));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Inspection of Alignment %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Connectivity map for skeleton drawing
% Let's take a look at the basic L2-norm fit
% line joining head to tip of left foot
Seg1    = [4, 3, 21, 2, 1, 13, 14, 15, 16];
% line from Spine base to right foot
Seg2    = [1, 17, 18, 19, 20];
% line from left hand to right hand
Seg3    = [22, 8, 7, 6, 5, 21, 9, 10, 11, 12, 24];
% left thumb Segment
Seg4    = [23, 7];
% right thumb Segment
Seg5    = [25, 11];

%%
% setup up first frame
K2      = reshape(K2Data(1,:),[3 25]);
SIMM    = reshape(SIMMData(1,:),[3 25]);

% estimate least squares solution
[s,R,T] = estimateSimTrans(K2,SIMM);

% transform Kinect2 data
K2      = s*R*K2 + repmat(T,[1 25]);

%Determine axis bounds by first taking the max/min of the max/min across 
%systems and then separating by X, Y, and Z.
tmpMax0 = max(SIMMData);
tmpMin0 = min(SIMMData);

tmpMax1 = max(K2Data);
tmpMin1 = min(K2Data);

tmpMax  = max([tmpMax0; tmpMax1]);
tmpMin  = max([tmpMin0; tmpMin1]);

% Separate by X, Y, and Z
Xmax    = max(tmpMax(1:3:end));
Xmin    = min(tmpMin(1:3:end));

Ymax    = max(tmpMax(2:3:end));
Ymin    = min(tmpMin(2:3:end));

Zmax    = max(tmpMax(3:3:end));
Zmin    = min(tmpMin(3:3:end));

FIG = figure('units','pixels','Position',[300 20 640 480],...
       'color','w',...
       'Resize','off');
   
set(FIG,'Renderer','zbuffer')

MainSubplot = subplot(1,3,2);
s0 = plot3(SIMM(1,:),SIMM(2,:),SIMM(3,:),'ok','clipping','off'); hold on
s1 = plot3(SIMM(1,Seg1),SIMM(2,Seg1),SIMM(3,Seg1),'k','clipping','off');
s2 = plot3(SIMM(1,Seg2),SIMM(2,Seg2),SIMM(3,Seg2),'k','clipping','off');
s3 = plot3(SIMM(1,Seg3),SIMM(2,Seg3),SIMM(3,Seg3),'k','clipping','off');
s4 = plot3(SIMM(1,Seg4),SIMM(2,Seg4),SIMM(3,Seg4),'k','clipping','off');
s5 = plot3(SIMM(1,Seg5),SIMM(2,Seg5),SIMM(3,Seg5),'k','clipping','off');

k0 = plot3(K2(1,:),K2(2,:),K2(3,:),'ob','clipping','off');
k1 = plot3(K2(1,Seg1),K2(2,Seg1),K2(3,Seg1),'b','clipping','off');
k2 = plot3(K2(1,Seg2),K2(2,Seg2),K2(3,Seg2),'b','clipping','off');
k3 = plot3(K2(1,Seg3),K2(2,Seg3),K2(3,Seg3),'b','clipping','off');
k4 = plot3(K2(1,Seg4),K2(2,Seg4),K2(3,Seg4),'b','clipping','off');
k5 = plot3(K2(1,Seg5),K2(2,Seg5),K2(3,Seg5),'b','clipping','off');

cameratoolbar
cameratoolbar('SetCoordSys','y')
set(MainSubplot,'xlim',[Xmin Xmax])
set(MainSubplot,'ylim',[Ymin Ymax])
set(MainSubplot,'zlim',[Zmin Zmax])
view(0,90)
axis equal
axis manual
axis off

SubplotSide = subplot(1,3,1);
copyobj(get(MainSubplot,'children'), SubplotSide);cameratoolbar
cameratoolbar('SetCoordSys','y')
set(SubplotSide,'xlim',[Xmin Xmax])
set(SubplotSide,'ylim',[Ymin Ymax])
set(SubplotSide,'zlim',[Zmin Zmax]), 
axis equal
axis manual
axis off

SubplotOH   = subplot(1,3,3);
copyobj(get(MainSubplot,'children'), SubplotOH);
cameratoolbar
cameratoolbar('SetCoordSys','y')
set(SubplotOH,'xlim',[Xmin Xmax])
set(SubplotOH,'ylim',[Ymin Ymax])
set(SubplotOH,'zlim',[Zmin Zmax]),
view(0,180)
axis equal
axis manual
axis off

TITLE = annotation('textbox',[.2 .8 .1 .1]);
    set(TITLE,'string',...
        sprintf('Aligning Kinect2 to SIMM: frame %d out of %d',1, numel(Frames)));
    set(TITLE,'FontSize',12)
    set(TITLE,'LineStyle','none')

% manually set some figure properties

%% Lock in properties and record
zoom off
% set(MainSubplot,'CameraPositionMode','manual')
% set(SubplotSide,'CameraPositionMode','manual')
% set(SubplotOH,'CameraPositionMode','manual')
frame            = getframe(FIG);

writer           = VideoWriter([VideoTitle '.mp4'],'MPEG-4');
writer.FrameRate = fps;
writer.Quality   = 50;
open(writer);
writeVideo(writer, frame);

for ii = 2:numel(Frames)
    K2      = reshape(K2Data(ii,:),[3 25]);
    SIMM    = reshape(SIMMData(ii,:),[3 25]);
    
    % estimate least squares solution
    [s,R,T] = estimateSimTrans(K2,SIMM);

    % transform Kinect2 data
    K2      = s*R*K2 + repmat(T,[1 25]);
    
    set(s0,'xdata',SIMM(1,:),'ydata',SIMM(2,:),'zdata',SIMM(3,:))
    set(s1,'xdata',SIMM(1,Seg1),'ydata',SIMM(2,Seg1),'zdata',SIMM(3,Seg1))
    set(s2,'xdata',SIMM(1,Seg2),'ydata',SIMM(2,Seg2),'zdata',SIMM(3,Seg2))
    set(s3,'xdata',SIMM(1,Seg3),'ydata',SIMM(2,Seg3),'zdata',SIMM(3,Seg3))
    set(s4,'xdata',SIMM(1,Seg4),'ydata',SIMM(2,Seg4),'zdata',SIMM(3,Seg4))
    set(s5,'xdata',SIMM(1,Seg5),'ydata',SIMM(2,Seg5),'zdata',SIMM(3,Seg5))
    axis([Xmin Xmax Ymin Ymax Zmin Zmax])
    axis equal

    set(k0,'xdata',K2(1,:),'ydata',K2(2,:),'zdata',K2(3,:))
    set(k1,'xdata',K2(1,Seg1),'ydata',K2(2,Seg1),'zdata',K2(3,Seg1))
    set(k2,'xdata',K2(1,Seg2),'ydata',K2(2,Seg2),'zdata',K2(3,Seg2))
    set(k3,'xdata',K2(1,Seg3),'ydata',K2(2,Seg3),'zdata',K2(3,Seg3))
    set(k4,'xdata',K2(1,Seg4),'ydata',K2(2,Seg4),'zdata',K2(3,Seg4))
    set(k5,'xdata',K2(1,Seg5),'ydata',K2(2,Seg5),'zdata',K2(3,Seg5))
    axis([Xmin Xmax Ymin Ymax Zmin Zmax])
    axis equal
    
    % update views
    copyobj(get(MainSubplot,'children'), SubplotSide)
    copyobj(get(MainSubplot,'children'),SubplotOH)
    
    set(TITLE, 'string', ...
        sprintf('Aligning Kinect2 to SIMM: frame %d out of %d',ii, numel(Frames)));
    drawnow
    
    frame            = getframe(FIG);
    writeVideo(writer, frame);
    
    cla(SubplotSide)
    cla(SubplotOH)
end

close(writer)