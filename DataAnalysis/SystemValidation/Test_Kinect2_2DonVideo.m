% TEST 2D camera projection data onto compressed video file
clear
clc

% test data folder
folder  = getenv('KinectData');
% test files
fileBase    = 'KinectTracking_2015-Apr-09-051107';
kinect2Data = [fileBase '_2D.bin'];
kinect2Vid  = [fileBase '.mp4'];

% Load Kinect2 data
fileInfo    = dir(sprintf('%s%s',folder,kinect2Data));
NJoints     = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(NJoints*2*SizeOfFloat);
fprintf(1,'Number of frames in Kinect2 file: %d.\n',NFrames);

% Compose Data as [2*NJoints NFrames] matrix and transpose
fid         = fopen(sprintf('%s%s',folder,kinect2Data),'r');
K2Data      = fread(fid,[2*NJoints NFrames],'float');
K2Data      = K2Data';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Inspection of 2D Data on Video %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% setup up first frame
K2        = reshape(K2Data(1,:),[2 25]);

T0 = tic;
fprintf(1,'Reading in videos frames...');
vidFrames = mmread(sprintf('%s%s',folder,kinect2Vid));
fprintf(1,'done in %9.6f seconds.\n',toc(T0));

tmp1 = flipud(squeeze(vidFrames.frames(1).cdata(:,:,1)));
tmp2 = flipud(squeeze(vidFrames.frames(1).cdata(:,:,2)));
tmp3 = flipud(squeeze(vidFrames.frames(1).cdata(:,:,3)));
tmp  = cat(3,tmp1,tmp2,tmp3);

V  = imshow(tmp); hold on

FIG = gcf;

k0 = plot(K2(1,:),K2(2,:),'ob','clipping','off');
k1 = plot(K2(1,Seg1),K2(2,Seg1),'b','clipping','off');
k2 = plot(K2(1,Seg2),K2(2,Seg2),'b','clipping','off');
k3 = plot(K2(1,Seg3),K2(2,Seg3),'b','clipping','off');
k4 = plot(K2(1,Seg4),K2(2,Seg4),'b','clipping','off');
k5 = plot(K2(1,Seg5),K2(2,Seg5),'b','clipping','off');

for ii = 2:NFrames
    K2       = reshape(K2Data(ii,:),[2 25]);
    
    tmp1 = flipud(squeeze(vidFrames.frames(ii).cdata(:,:,1)));
    tmp2 = flipud(squeeze(vidFrames.frames(ii).cdata(:,:,2)));
    tmp3 = flipud(squeeze(vidFrames.frames(ii).cdata(:,:,3)));
    tmp  = cat(3,tmp1,tmp2,tmp3);
    
    set(V,'CData',tmp);
    
    set(k0,'xdata',K2(1,:),'ydata',K2(2,:))
    set(k1,'xdata',K2(1,Seg1),'ydata',K2(2,Seg1))
    set(k2,'xdata',K2(1,Seg2),'ydata',K2(2,Seg2))
    set(k3,'xdata',K2(1,Seg3),'ydata',K2(2,Seg3))
    set(k4,'xdata',K2(1,Seg4),'ydata',K2(2,Seg4))
    set(k5,'xdata',K2(1,Seg5),'ydata',K2(2,Seg5))
    
    drawnow
end