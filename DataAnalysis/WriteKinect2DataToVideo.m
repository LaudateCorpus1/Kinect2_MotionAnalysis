function WriteKinect2DataToVideo(kinect2BinFile)
%WRITEKINECT2DATATOVIDEO : 
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
% a mp4 video file in the same location and named according to the kinect2
% binary file.
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
[filepath, filename, ext] = fileparts(kinect2BinFile);

% check for underscores in filename
tmp = regexp(filename,'[_]');
if ~isempty(tmp)
    strName = filename;
    strName(tmp) = ' ';
else
    strName = filename;
end

if ~strcmp(ext, '.bin')
    error('%s: The file ''%s'' is not a binary file.', mfname, kinect2BinFile);
end

if isempty(filepath)
    filepath = pwd;
end

% Load Kinect2 data
fid         = fopen(kinect2BinFile,'r');
NJoints     = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(NJoints*3*SizeOfFloat);

if (NFrames - floor(NFrames)) > 0
    error(['%s: Incorrect File: The computed number of frames in ''%s'' is '...
           'not an integer value.'], mfname, kinect2BinFile);
end
fprintf(1,'Number of frames in Kinect2 file %s: %d.\n',kinect2BinFile, NFrames);

% Compose Data as [3*NJoints NFrames] matrix and transpose
K2Data      = fread(fid,[3*NJoints NFrames],'float');
K2Data      = K2Data';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Connectivity map for skeleton drawing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
K2      = reshape(K2Data(1,:),[3 NJoints]);

% Determine axis bounds
Maxs   = max(K2,[],2);
Mins   = min(K2,[],2);
factor = 1.10;
% Separate by X, Y, and Z
Xmax = factor*Maxs(1);
Xmin = factor*Mins(1);

Ymax = factor*Maxs(2);
Ymin = factor*Mins(2);

Zmax = factor*Maxs(3);
Zmin = factor*Mins(3);

FIG = figure('units','pixels','Position',[300 20 640 480],'color','w');

k0 = plot3(K2(1,:),K2(2,:),K2(3,:),'ok','markerfacecolor','k'); hold on
k1 = plot3(K2(1,Seg1),K2(2,Seg1),K2(3,Seg1),'b');
k2 = plot3(K2(1,Seg2),K2(2,Seg2),K2(3,Seg2),'b');
k3 = plot3(K2(1,Seg3),K2(2,Seg3),K2(3,Seg3),'b');
k4 = plot3(K2(1,Seg4),K2(2,Seg4),K2(3,Seg4),'b');
k5 = plot3(K2(1,Seg5),K2(2,Seg5),K2(3,Seg5),'b');

cameratoolbar
cameratoolbar('SetCoordSys','y')
axis([Xmin Xmax Ymin Ymax Zmin Zmax])
axis equal
axis off

camtarget(K2(:,1)');

TITLE = annotation('textbox',[.6 .1 .1 .1]);
    set(TITLE,'string',...
        sprintf('%s:\nframe %6d out of %6d',strName, 1, NFrames));
    set(TITLE,'FontSize',12)
    set(TITLE,'LineStyle','none')

    
fprintf(1,'Set final camera angle prior to writing video...\n');
keyboard

% Lock in properties and record
frame            = getframe(FIG);
VideoTitle       = filename;
writer           = VideoWriter([VideoTitle '.mp4'],'MPEG-4');
writer.FrameRate = 30;
writer.Quality   = 50;
open(writer);
writeVideo(writer, frame);

for ii = 2:NFrames
    K2      = reshape(K2Data(ii,:),[3 NJoints]);
    
    set(k0,'xdata',K2(1,:),'ydata',K2(2,:),'zdata',K2(3,:))
    set(k1,'xdata',K2(1,Seg1),'ydata',K2(2,Seg1),'zdata',K2(3,Seg1))
    set(k2,'xdata',K2(1,Seg2),'ydata',K2(2,Seg2),'zdata',K2(3,Seg2))
    set(k3,'xdata',K2(1,Seg3),'ydata',K2(2,Seg3),'zdata',K2(3,Seg3))
    set(k4,'xdata',K2(1,Seg4),'ydata',K2(2,Seg4),'zdata',K2(3,Seg4))
    set(k5,'xdata',K2(1,Seg5),'ydata',K2(2,Seg5),'zdata',K2(3,Seg5))
    
    % Determine axis bounds
    Maxs   = max(K2,[],2);
    Mins   = min(K2,[],2);
    % Separate by X, Y, and Z
    Xmax = factor*Maxs(1);
    Xmin = factor*Mins(1);

    Ymax = factor*Maxs(2);
    Ymin = factor*Mins(2);

    Zmax = factor*Maxs(3);
    Zmin = factor*Mins(3);
    axis([Xmin Xmax Ymin Ymax Zmin Zmax])
    
    camtarget(K2(:,1)');
    
    axis equal
    axis auto
    
    set(TITLE, 'string', ...
        sprintf('%s:\nframe %6d out of %6d',strName, ii, NFrames));
    drawnow
    
    frame = getframe(FIG);
    writeVideo(writer, frame);
    
end

close(writer)
close(FIG)
