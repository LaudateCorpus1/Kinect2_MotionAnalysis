function kinect2QualityControl(FileName)
% kinect2QualityControl: This code takes in a file name formatted according
% to the Kinect2Acq system and provides a windows media player instance for
% conveniently assessing the quality of the data
%
%%% INPUTS %%%
%
% FileName: a string path to the data (no extension
%           -if none supplied, a file system interface is opened  
%
% author: John D. Long II, PhD, contact: jlong29@gmail.com

mname = mfilename;

if ~ispc
    error('%s::UseError::This code requires ActiveX and is therefore Windows based',mname);
end

if nargin < 1 || isempty(FileName)
    [vidFile pathName]  = uigetfile('*.*','Please select a video to validate'); %#ok<NCOMMA>
    [~,fileName]        = fileparts(vidFile);
    kinect2File         = [fileName '_3D.bin'];  
else
    % check file path
    chkFile                     = dir([FileName '*']);
    if isempty(chkFile)
        error('%s::InputError::FileNotFound',mname)
    else
        % Check if all necessary files are available
        fileChks = {
                    [FileName '.mp4'],...
                    [FileName '_3D.bin']
                    };
        if sum(ismember({chkFile(:).name},fileChks)) < numel(fileChks)
            error('%s::InputError::Cannot find .mp4 and/or _3d.bin file');
        else
            vidFile     = fileChks{1};
            kinect2File = fileChks{2};
            
            fullPath    = which(vidFile);
            pathName    = [fileparts(fullPath) filesep];
        end
    end
end

% Access monitor parameters and initialize figure layout
set(0,'Units','pixels')
scrn        = get(0,'MonitorPositions');
scrnsz      = scrn(1,:);
aspMonitor  = scrnsz(3)/scrnsz(4);

if isequal(aspMonitor, 16/9)   % 1080 aspect ratio
    pos = scrnsz(1,:);
elseif aspMonitor > 16/9
    pos = [scrnsz(1:3) floor(scrnsz(3)*(9/16))];
else
    pos = scrn(1,:);
end

%%%%%%%%%%%%%%%%%
%%% LOAD DATA %%%
%%%%%%%%%%%%%%%%%

% set up video file placement in figure in pixels coords (activex requires)
vidInfo   = VideoReader(fullfile(pathName, vidFile));
fdur      = 1/vidInfo.FrameRate;

% Load Kinect2 data
fileInfo    = dir(fullfile(pathName, kinect2File));
fid         = fopen(fullfile(pathName, kinect2File),'r');
NJoints     = 25;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(NJoints*3*SizeOfFloat);

if (NFrames - floor(NFrames)) > 0
    error(['%s: Incorrect File: The computed number of frames in ''%s'' is '...
           'not an integer value.'], mfname, kinect2File);
end
fprintf(1,'Number of frames in Kinect2 file %s: %d.\n',kinect2File, NFrames);

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

% Test markers for labeling nodes
Test = [4, 3, 21, 2, 1, 5, 6, 7];

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

% Parameters for Display
aspVid    = vidInfo.Height/vidInfo.Width;

vidPosX   = ceil(0.05*pos(3));
vidPosY   = ceil(0.2*pos(4));
vidWidth  = pos(3)/2;
vidHeight = aspVid*vidWidth;

% Initialize Figure
H1       = figure('units','pixels','Position',pos,'color','w');
Actx     = actxcontrol('WMPlayer.OCX.7', [vidPosX vidPosY vidWidth vidHeight], H1);
Actx.URL = [pathName vidFile];
waitfor(Actx,'playState');
while strcmp(Actx.playState,'wmppsPlaying')
    disp(Actx.playState)
    Actx.controls.pause;
    set(Actx.controls,'currentPosition',0);
end

set(H1,'units','normalized');

% set uicontrol button controls
RestartBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.2 0.11 0.05 0.05], ...
                         'String', 'o', ...
                         'Callback', @(src,event)RestartBtnCallback());


PlayBtn    = uicontrol(H1,'units','normalized', ...
                         'Position', [0.3 0.11 0.05 0.05], ...
                         'String', 'Play', ...
                         'Callback', @(src,event)PlayBtnCallback());

PauseBtn   = uicontrol(H1,'units','normalized', ...
                         'Position', [0.25 0.11 0.05 0.05], ...
                         'String', 'Pause', ...
                         'Callback', @(src,event)PauseBtnCallback());

StepFwdBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.35 0.11 0.05 0.05], ...
                         'String', '>', ...
                         'Callback', @(src,event)StepFwdBtnCallback());

% set keyPress function for hotkeys
set(H1,'KeyPressFcn',@action_hotkeys);
set(H1,'CloseRequestFcn',@closeFunction);
% Set and save guidata
Gui_Data.Actx = Actx;
guidata(H1,Gui_Data);

% set kinematic data frame
k2PosX   = ceil(0.55*pos(3));
k2PosY   = ceil(0.2*pos(4));
k2Width  = pos(3)/2;
k2Height = aspVid*vidWidth;

set(H1,'units','pixels');
ax1 = axes('units','pixels','Position',[k2PosX k2PosY k2Width k2Height],...
            'color','w','parent',H1);

k0 = plot3(K2(1,:),K2(2,:),K2(3,:),'ok','markerfacecolor','k'); hold on
k1 = plot3(K2(1,Seg1),K2(2,Seg1),K2(3,Seg1),'b');
k2 = plot3(K2(1,Seg2),K2(2,Seg2),K2(3,Seg2),'b');
k3 = plot3(K2(1,Seg3),K2(2,Seg3),K2(3,Seg3),'b');
k4 = plot3(K2(1,Seg4),K2(2,Seg4),K2(3,Seg4),'b');
k5 = plot3(K2(1,Seg5),K2(2,Seg5),K2(3,Seg5),'b');
k6   = plot3(K2(1,Test),K2(2,Test),K2(3,Test),'or');

cameratoolbar
cameratoolbar('SetCoordSys','y')
axis([Xmin Xmax Ymin Ymax Zmin Zmax])
axis equal
axis off

camtarget(K2(:,1)');

strName                      = kinect2File;
strName(regexp(strName,'_')) = ' ';
TITLE = annotation('textbox',[.6 .8 .1 .1]);
    set(TITLE,'string',...
        sprintf('%s:\nFrame %4d of %4d',strName, 1, NFrames));
    set(TITLE,'FontSize',12)
    set(TITLE,'LineStyle','none') 

FrameCounter = 1;
Replay       = true;
while Replay
    while FrameCounter <= NFrames
        T0 = tic;
        
        if ~isvalid(H1)
            return;
        end
        
        % update video
        StepFwdBtnCallback();
        
        % update kinect data
        K2      = reshape(K2Data(FrameCounter,:),[3 NJoints]);

        set(k0,'xdata',K2(1,:),'ydata',K2(2,:),'zdata',K2(3,:))
        set(k1,'xdata',K2(1,Seg1),'ydata',K2(2,Seg1),'zdata',K2(3,Seg1))
        set(k2,'xdata',K2(1,Seg2),'ydata',K2(2,Seg2),'zdata',K2(3,Seg2))
        set(k3,'xdata',K2(1,Seg3),'ydata',K2(2,Seg3),'zdata',K2(3,Seg3))
        set(k4,'xdata',K2(1,Seg4),'ydata',K2(2,Seg4),'zdata',K2(3,Seg4))
        set(k5,'xdata',K2(1,Seg5),'ydata',K2(2,Seg5),'zdata',K2(3,Seg5))
        
        set(k6,'xdata',K2(1,Test),'ydata',K2(2,Test),'zdata',K2(3,Test))
        
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
            sprintf('%s:\nFrame %4d of %4d',strName, FrameCounter, NFrames));
        drawnow
        pause(fdur-toc(T0))

        % Update FrameCounter
        FrameCounter = FrameCounter + 1;
    end
    RestartBtnCallback();
end

    %%%%%%%%%%%%%%%%%
    %%% CALLBACKS %%%
    %%%%%%%%%%%%%%%%%
    % Restart Video button
    function RestartBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.play;GUI_Data.Actx.controls.pause
        set(GUI_Data.Actx.controls,'currentPosition',0);
        FrameCounter = 1;
    end
    
    % play button
    function PlayBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.play;
    end
    
    % pause button
    function PauseBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.pause;
    end
    
    % step forward button
    function StepFwdBtnCallback()
        if isvalid(H1)
            GUI_Data = guidata(H1);
            GUI_Data.Actx.controls.pause;
            GUI_Data.Actx.controls.step(1);
            disp(GUI_Data.Actx.controls.currentPosition)
            disp(GUI_Data.Actx.playState)
        end
    end
    
    % HotKeys for speeding up workflow
    function action_hotkeys(~,evnt)
        % Select Actions
        if length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'z')
            
            RestartBckBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'x')
            
            PauseBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'c')
            
            PlayBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'v')
            
            StepFwdBtnCallback()
            
        end
    end
    
    % figure close function
    function closeFunction(src,callbackdata) %#ok<INUSD>
        % clean up
        GUI_Data = guidata(H1);
        GUI_Data.Actx.close;
        delete(GUI_Data.Actx);
        delete(H1);
        return;
    end
end