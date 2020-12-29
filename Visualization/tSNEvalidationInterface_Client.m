function tSNEvalidationInterface_Client(videoFile,tSNEmap, Embedding)

if nargin < 1 || isempty(videoFile)
    [videoFile pathName] = uigetfile('*.*','Please select a video to validate'); %#ok<NCOMMA>
    videoFile            = [pathName videoFile];
end

% check file path
[pathName,fileName,ext] = fileparts(videoFile);
if isempty(pathName)
    chk = which(videoFile);
    [pathName,fileName,ext] = fileparts(chk);
else
    chk = dir([pathName filesep fileName ext]);
    chk = chk.name;    
end

if isempty(chk)
    error('Input Error: the videoFile entered is either incomplete or not in the MATLAB path.')
else
    fprintf(1,'Loading Video File %s...\n',videoFile);
    pathName  = [pathName filesep];
    videoFile = [fileName ext];
end

% Access monitor parameters and initialize figure layout
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

% set up video file placement in figure in pixels coords (activex requires)
vidInfo = VideoReader([pathName videoFile]);
aspVid  = vidInfo.Width/vidInfo.Height;
clear vidInfo

vidPosX   = ceil(0.5*scrnsz(3));
vidPosY   = ceil(scrnsz(4)/4);
vidWidth  = ceil(scrnsz(3)/2);
vidHeight = floor(vidWidth/aspVid);

% Initialize Figure
H1 = figure('units','pixels','Position',pos,'color','w');

% Activex interface with Windows Media Player (TO DO: linux, vlc)
Actx     = actxcontrol('WMPlayer.OCX.7', [vidPosX vidPosY vidWidth vidHeight], H1);
Actx.URL = [pathName videoFile];
waitfor(Actx,'playState');
Actx.controls.pause;

% switch to normalized coordinates for setting up the rest of the figure
set(H1,'units','normalized');

% set uicontrol button controls
PlayBtn    = uicontrol(H1,'units','normalized', ...
                         'Position', [0.75 0.11 0.05 0.05], ...
                         'String', 'Play', ...
                         'Callback', @(src,event)PlayBtnCallback());

PauseBtn   = uicontrol(H1,'units','normalized', ...
                         'Position', [0.70 0.11 0.05 0.05], ...
                         'String', 'Pause', ...
                         'Callback', @(src,event)PauseBtnCallback());

StepFwdBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.80 0.11 0.05 0.05], ...
                         'String', '>', ...
                         'Callback', @(src,event)StepFwdBtnCallback());

StepBckBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.65 0.11 0.05 0.05], ...
                         'String', '<', ...
                         'Callback', @(src,event)StepBckBtnCallback());

FastFwdBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.85 0.11 0.05 0.05], ...
                         'String', '>>', ...
                         'Callback', @(src,event)FastFwdBtnCallback());

                     
FastBckBtn = uicontrol(H1,'units','normalized', ...
                         'Position', [0.60 0.11 0.05 0.05], ...
                         'String', '<<', ...
                         'Callback', @(src,event)FastBckBtnCallback());

% set keyPress function for hotkeys
set(H1,'KeyPressFcn',@action_hotkeys);
set(H1,'CloseRequestFcn',@closeFunction);
% Set and save guidata
Gui_Data.Actx = Actx;
guidata(H1,Gui_Data);

% Actx.controls.step(1)          % steps forward 1 frame
% Actx.controls.step(-1)         % steps back 1 second
% Actx.controls.currentPosition; % can be mapped to frames
    
    %%%%%%%%%%%%%%%%%
    %%% CALLBACKS %%%
    %%%%%%%%%%%%%%%%%
    
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
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.pause;
        GUI_Data.Actx.controls.step(1);
    end
    
    % step backward button
    function StepBckBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.pause;
        GUI_Data.Actx.controls.step(-1);
    end
    
    % fast forward button
    function FastFwdBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.pause;
        GUI_Data.Actx.controls.fastForward;
    end
    
    % fast reverse button
    function FastBckBtnCallback()
        GUI_Data = guidata(H1);
        GUI_Data.Actx.controls.pause;
        GUI_Data.Actx.controls.fastReverse;
    end

    % HotKeys for speeding up workflow
    function action_hotkeys(~,evnt)
        % Select Actions
        if length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'z')
            
            StepBckBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'x')
            
            PauseBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'c')
            
            PlayBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'v')
            
            StepFwdBtnCallback()
            
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'d')
            % Move node
            disp('Action: Move node')
            set(tl2,'string','Move Node','fontsize',12)
            action = 'drag';
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'r')
            % Rotate data_frame
            disp('Action: Rotate data_frame')
            set(tl1,'string',[{'Rotate Data Frame:'};{'Click in 3d plot to begin.'}])
            action = 'rotate';
        elseif length(evnt.Modifier) == 1 && strcmp(evnt.Modifier{:},'control') && strcmp(evnt.Key,'u')
            % Declare up vector
            disp('Action: Declare Up Vector')
            set(tl2,'string','Define Global Orientation Axis')
            action = 'up_vec';
            rcentroid = rot*centroid;
            set(cp2,'xdata',rcentroid(1),'ydata',rcentroid(2),'visible','on')
            set(a2,'visible','on')
            set(up_ptr,'visible','on')
        end
    end

    % figure close function
    function closeFunction(src,callbackdata) %#ok<INUSD>
        % clean up
        GUI_Data = guidata(H1);
        GUI_Data.Actx.close;
        delete(GUI_Data.Actx);
        delete(H1);
    end
end