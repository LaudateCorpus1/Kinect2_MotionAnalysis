%%% Data Analysis of 10.10.2014 data %%%
clc
clear

% folder containing the data
DataDir1 = [getenv('Kinect2') filesep 'Data\Data10.10.2014'];
DataDir2 = [getenv('Kinect2') filesep 'Data\Data06.14.2015'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set which variables to include in the analysis
Vars          = 1:75;

% DOF eliminated by normalization: base of spine, z-coords of HipRight and
% SpineMid
exclude       = [1:3 6 51];
Vars(exclude) = [];

% lowpass filters kinematics
lowpassK   = false;
Fs         = 30;   % sampling rate
Lp         = 5;    % low pass cutoff

% PCA denoising
PC_denoise = true;
InspectPC  = true;  % inspect PC component threshold
NoiseChk1  = false; % diagnostic figure for examining noise in PC time-series

% Wavelet decomposition
fMin          = 0.5;
fMax          = Lp;
nbins         = 25;
InspectW      = false;   % inspect wavelets
saveGroupData = true;    % save GroupData data to file
NoiseChk2     = false;   % diagnostic figure of examining noise in wavelet spectra

% set parameters for converting distance matrix to probability matrix
u            = 64;       % transition prob target value ~Ndatapoints/Nbehaviors
u_training   = u;        % perplexity on training set
tol          = 1.0e-5;   % tolerance for target value approximation

% tSNE
Nsub              = 20000;   % the number of samples in use when running tSNE
stepSz            = 1;      % step size for temporal subsampling
tnumPoints        = 20000;   % samples per session to use in training set
minTemplateLength = 10;      % minimum watershed region size

% Joint Type Enumeration (need to store as a file)
keys   = 1:25;
values = {...
            'SpineBase';
            'SpineMid';
            'Neck';
            'Head';
            'ShoulderLeft';
            'ElbowLeft';
            'WristLeft';
            'HandLeft';
            'ShoulderRight';
            'ElbowRight';
            'WristRight';
            'HandRight';
            'HipLeft';
            'KneeLeft';
            'AnkleLeft';
            'FootLeft';
            'HipRight';
            'KneeRight';
            'AnkleRight';
            'FootRight';
            'SpineShoulder';
            'HandTipLeft';
            'ThumbLeft';
            'HandTipRight';
            'ThumbRight'};

JointMap = containers.Map(keys,values);

%% 0) Clean and Compose data

% For the moment, I'm assuming all files have been normalized and converted
% to .mat files. This can be generalized.

% Load normalized .mat files containing kinect2 data and compose 
% into a structure with fields:
%   .Data: a [NFrames x Nnodes] matrix of kinematic data
%   .NFrames: the number of frames in the kinect2 data file
%   .fileInfo: a string detailing which file was loaded

% compose file list
Files = ...
    {[DataDir1 '\Anna\Sitting_Kinect_normalized.mat'];
     [DataDir1 '\Anna\Standing_Kinect_normalized.mat'];
     [DataDir1 '\Cory\Sitting_Kinect_normalized.mat'];
     [DataDir1 '\Cory\Standing_Kinect_normalized.mat'];
     [DataDir1 '\Kathleen\Sitting_Kinect_normalized.mat'];
     [DataDir1 '\Kathleen\Standing_Kinect_normalized.mat'];
     [DataDir1 '\Victor\Sitting_Kinect_normalized.mat'];
     [DataDir1 '\Victor\Standing_Kinect_normalized.mat']};

Nepochs = numel(Files);

% Compose structure layout for group data
GroupData = struct('sessions',cell(1,1),'params',cell(1,1));
tmpStruct = struct('Trials',   cell(Nepochs, 1));

GroupData.params.fileInfo = Files;

% temporary structure of composing file data
K2Data  = struct('Data',     cell(Nepochs,1), ...
                 'wavelets', cell(Nepochs,1), ...
                 'NFrames',  cell(Nepochs,1), ...
                 'fileInfo', cell(Nepochs,1));

for ii = 1:Nepochs
    
    % attempt to load data
    try
        tmp        = load(Files{ii});
    catch DamnErr
        fprintf(1,'\t\tError: %s\n',DamnErr.message);
        return;
    end
    
    % transfer data structure
    K2Data(ii).Data     = tmp.K2Data.Data;
    K2Data(ii).NFrames  = tmp.K2Data.NFrames;
    K2Data(ii).fileInfo = tmp.K2Data.fileInfo;
    
    % Transpose data to [NFrames x Nnodes]
    K2Data(ii).Data = K2Data(ii).Data';
    
    % exclude variables
    K2Data(ii).Data = K2Data(ii).Data(:,Vars);
end

%% 1) Lowpass filter and run PCA on complete data set

% log parameters
GroupData.params.general.Fs   = Fs;
GroupData.params.general.Vars = Vars;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRE-PROCESS DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% low pass filters kinematic data
if lowpassK
    fprintf(1,'Low pass filtering kinematic data...\n');
    
    % log parameters
    GroupData.params.lowPass.Lp = Lp;
    fMax                        = Lp;
    
    % build filter
    d  = fdesign.lowpass('Fp,Fst,Ap,Ast',Lp,Lp+2,1,80,Fs);
    Hd = design(d,'butter');
    
    GroupData.params.lowPass.sosMatrix   = Hd.sosMatrix;
    GroupData.params.lowPass.Scalevalues = Hd.ScaleValues;
    
    for ii = 1:numel(K2Data)
        fprintf(1,'File: %s\n',K2Data(ii).fileInfo);
        tmp = K2Data(ii).Data;
        tmp = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmp);
        K2Data(ii).Data = tmp;
    end
end

% Denoise using PCA
if PC_denoise
    fprintf(1,'Running PCA on the data...');
    
    [K2Data, PCs, lastEigen] = Kinect2PCs(K2Data,InspectPC);
    
    % log parameters
    GroupData.params.PCA.lastEigen = lastEigen;
    GroupData.params.PCA.meanData  = PCs.meanData;
    GroupData.params.PCA.stdData   = PCs.stdData;
    GroupData.params.PCA.projM     = PCs.projM;
    GroupData.params.PCA.eigV      = PCs.eigV;
    fprintf(1,'done.\n');
else
    % working from raw data
    fprintf(1,'PCA not run...\n');
    % log parameters
    GroupData.params.PCA.lastEigen = size(K2Data(1).Data,2);
    GroupData.params.PCA.meanData  = nan;
    GroupData.params.PCA.stdData   = nan;
    GroupData.params.PCA.projM     = nan;
    GroupData.params.PCA.eigV      = nan;
end

% RECONSTRUCTION EXAMPLE
% tmp0 = K2Data(1).Data(1,:);
% chk0 = reshape(tmp0,[3 25]);
%
% chk1 = [nan(3,1);sum(repmat(tmp1(1:40),[72 1]).*PCs.projM(:,1:40),2)]'.*PCs.stdData+PCs.meanData
% chk1 = reshape(chk1,[3 25]);
%
% figure, plot3(chk0(1,:),chk0(2,:),chk0(3,:),'.k')
% hold on, plot3(chk1(1,:),chk1(2,:),chk1(3,:),'oc')
% axis equal,cameratoolbar, grid on

%% NoiseChk1 : diagnostic figure for examining noise in PC time-series
if NoiseChk1
    % Useful diagnostic figure for determining noise ceiling of raw data
    if exist('data','var')
        if size(data,1) ~= sizeChk0
            fprintf(1,'Loading flat data matrix...\n');
            data = cat(1,K2Data(:).Data);
            sizeChk0 = size(data,1);
        else
            fprintf('Data matrix already loaded...\n');
        end
    else
        fprintf(1,'Loading flat data matrix...\n');
        data = cat(1,K2Data(:).Data);
        sizeChk0 = size(data,1);
    end
    
    % calculate axis limits
    maxPC  = max(data(:));
    minPC  = min(data(:));
    offset = 1;
    winSz  = 500;
    
    % Inspect the PCs up to one passed the last significant PC
    FIGCHK = figure('units','pixels','position',[400 150 1000 500],'color','w');
    plot(data(:,GroupData.params.PCA.lastEigen),'k')
    hold on,
    plot(data(:,1:GroupData.params.PCA.lastEigen))
    z = zoom;
    setAxesZoomMotion(z,gca,'horizontal');
    p = pan;
    setAxesPanMotion(p,gca,'horizontal');
    set(gca,'ylim',[minPC maxPC])
    set(gca,'xlim',[offset offset+winSz])
    ylabel('(a.u.)')
    xlabel('Time (samples)')
    title('Time-series Data: PC or Raw ')
    
end

%% 2) Apply wavelet decomposition to kinect2 data

% OUTPUT: Decompose into wavelets and save into compute directory
waveletFeatures = Kinect2Wavelets(K2Data,'fMin',fMin,'fMax',fMax,'nbins',nbins,'Inspect',InspectW);

% scale wavelets to mad = 1 and min per wavelet = 1e-5
    
% compute median statistics for this session
tmp     = cat(1,waveletFeatures(:).wavelets.power);
medData = nanmedian(tmp);
madData = mad(tmp);

for jj = 1:numel(waveletFeatures.wavelets)
    data       = waveletFeatures.wavelets(jj).power;

    data       = (data - repmat(medData,[size(data,1) 1]))./...
                    repmat(1.4826*madData,[size(data,1) 1]);
    data       = data - repmat(min(data)-1e-5,[size(data,1) 1]);

    waveletFeatures.wavelets(jj).power = data;
end
    
GroupData.params.wavelets = waveletFeatures.params;

% compose group data: TO DO partition by sessions
for ii = 1:Nepochs
    K2Data(ii).wavelets = waveletFeatures.wavelets(ii);
end
clear waveletFeatures

GroupData.sessions = K2Data;
clear K2Data

% log tSNE parameters
GroupData.params.tSNE.Nsub     = Nsub;
GroupData.params.tSNE.u        = u;
GroupData.params.tSNE.tol      = tol;

% compose Berman parameter structure
% set up parameters structure needed by embedding routine
parameters.trainingSetSize     = Nsub;
parameters.training_numPoints  = tnumPoints;
parameters.numProcessors       = 6;
parameters.pcaModes            = GroupData.params.PCA.lastEigen;
parameters.numPeriods          = nbins;
parameters.samplingFreq        = GroupData.params.general.Fs;
parameters.minF                = GroupData.params.wavelets.fMin;
parameters.maxF                = GroupData.params.wavelets.fMax;
parameters.omega0              = GroupData.params.wavelets.omega;
parameters.perplexity          = u;
parameters.training_perplexity = u_training;
parameters.sigmaTolerance      = tol;
parameters.maxNeighbors        = max([4*u 200]);    % don't let perplexity > maxNeighbors
parameters.numProjections      = GroupData.params.PCA.lastEigen*nbins;
parameters.kdNeighbors         = max([log2(u) 10]);
parameters.minTemplateLength   = minTemplateLength;

parameters.basisImage          = nan;

parameters                     = setRunParameters(parameters);

parameters.stepSz              = stepSz;

computeFolder = [getenv('Kinect2') filesep 'Data\Data10.10.2014' filesep 'DataProcessing'];

if saveGroupData
    fprintf(1,'Saving GroupData to %s...',computeFolder);
    save([computeFolder filesep 'GroupData_test.mat'],'GroupData')
    fprintf(1,'done.\n');
    
    fprintf(1,'Saving parameters for run to %s...',computeFolder);
    save([computeFolder filesep 'parameters_test.mat'],...
            'parameters')
    fprintf(1,'done.\n');
end

%% NoiseChk2 : diagnostic figure of examining noise in wavelet spectra
% Inspect max power at each frequency across all time points to find
% separation between behavioral and process (aka sensor and model fitting)
% dyanmics.

clear data
if NoiseChk2
    
    if exist('Pow','var')
        if size(Pow,1) ~= sizeChk0
            fprintf(1,'Loading flat Pow matrix...\n');
            tmp = cat(1,GroupData.sessions(:).wavelets);
            Pow      = cat(1,tmp(:).power);
            sizeChk0 = size(Pow,1);
            
            % access parameters
            nbins   = GroupData.params.wavelets.nbins;
            
            % for excluding the highest freqenc(ies)
            offset = 1;
            
            % first pass to set  bounds
            MAXVAL0 = zeros(size(Pow,2)/nbins,1);
            MAXVAL1 = zeros(size(Pow,2)/nbins,1);
            
            for ii = 1:size(Pow,2)/nbins
                chk  = Pow(:,nbins*(ii-1)+1:nbins*ii-offset);

                chk0 = max(chk);
                chk0 = chk0/min(chk0);

                MAXVAL0(ii) = max(chk0);

                chk1 = median(chk);
                chk1 = chk1/min(chk1);

                MAXVAL1(ii) = max(chk1);
            end
            MAXVAL0 = max(MAXVAL0);
            MINVAL0 = 0;

            MAXVAL1 = max(MAXVAL1);
            MINVAL1 = 0;
            
        else
            fprintf('Data matrix already loaded...\n');
        end
    else
        fprintf(1,'Loading flat Pow matrix...\n');
        tmp = cat(1,GroupData.sessions(:).wavelets);
        Pow      = cat(1,tmp(:).power);
        sizeChk0 = size(Pow,1);

        % access parameters
        nbins   = GroupData.params.wavelets.nbins;

        % for excluding the highest freqenc(ies)
        offset = 1;

        % first pass to set  bounds
        MAXVAL0 = zeros(size(Pow,2)/nbins,1);
        MAXVAL1 = zeros(size(Pow,2)/nbins,1);

        for ii = 1:size(Pow,2)/nbins
            chk  = Pow(:,nbins*(ii-1)+1:nbins*ii-offset);

            chk0 = max(chk);
            chk0 = chk0/min(chk0);

            MAXVAL0(ii) = max(chk0);

            chk1 = median(chk);
            chk1 = chk1/min(chk1);

            MAXVAL1(ii) = max(chk1);
        end
        MAXVAL0 = max(MAXVAL0);
        MINVAL0 = 0;

            MAXVAL1 = max(MAXVAL1);
            MINVAL1 = 0;
            
    end

    FIGCHK = figure('units','pixels','position',[400 150 1000 500],'color','w');
    for ii = 1:size(Pow,2)/nbins
        Node  = JointMap(ceil(Vars(ii)/3));
        chk   = mod(Vars(ii),3);
        if chk == 1
            Coord = 'X';
        elseif chk == 2
            Coord = 'Y';
        else
            Coord = 'Z';
        end
        
        chk    = Pow(:,nbins*(ii-1)+1:nbins*ii-offset);
        
        tmp0 = max(chk);
        tmp0 = tmp0/min(tmp0);
        
        s1 = subplot(1,2,1); 
        semilogy(GroupData.params.wavelets.freqs(1:end-offset),tmp0,'b','linewidth',1.5)
        title(sprintf('Maximum power across time points for Variable %d: %s, %s coord',ii, Node, Coord)),
        ylabel('Power Ratio (a.u.)')
        xlabel('Freqs (cycles/second')
        
        tmp1 = median(chk);
        tmp1 = tmp1/min(tmp1);
        
        s2 = subplot(1,2,2);
        semilogy(GroupData.params.wavelets.freqs(1:end-offset),tmp1,'k','linewidth',1.5)
        title(sprintf('Median power across time points for Variable %d: %s, %s coord',ii, Node, Coord)),
        ylabel('Power Ratio (a.u.)')
        xlabel('Freqs (cycles/second')
        
        set(s1,'ylim',[1 MAXVAL0])
        set(s1,'xlim',[GroupData.params.wavelets.freqs(1) GroupData.params.wavelets.freqs(end-offset)])
        
        set(s2,'ylim',[1 MAXVAL1])
        set(s2,'xlim',[GroupData.params.wavelets.freqs(1) GroupData.params.wavelets.freqs(end-offset)])
        
        pause
    end
    close(FIGCHK);
end

%% 3) Compose tSNE sample set across Files

fprintf(1,'Running ComposeTSNEsamples.m...\n');
try
    keep GroupData parameters
catch
    clear
end
curDir = pwd;

DEBUG_SAMPLES = 0;

% directory for caching variables
computeFolder = [getenv('Kinect2') filesep 'Data\Data10.10.2014' filesep 'DataProcessing'];

if ~exist('GroupData','var')
    fprintf(1,'Attempting to load GroupData...\n');
    
    cd(computeFolder)
    try
        TestFile = 'GroupData_test.mat';
        tmp      = dir(TestFile);
        fprintf(1,'Group Data Wavelets created at %s.\n',tmp.date);
        fprintf(1,'Importing wavelet features...');
        load([computeFolder filesep TestFile],'GroupData')
        fprintf(1,'done.\n');
        cd(curDir)
    catch DamnErr
        cd(curDir)
        error('Error: %s\n',DamnErr.message);
    end
end

if ~exist('parameters','var')
    fprintf(1,'Attempting to load parameters...\n');
    
    cd(computeFolder)
    try
        TestFile = 'parameters_test.mat';
        tmp      = dir(TestFile);
        fprintf(1,'Group Data parameters created at %s.\n',tmp.date);
        fprintf(1,'Importing parameters features...');
        load([computeFolder filesep TestFile],'parameters')
        fprintf(1,'done.\n');
        cd(curDir)
    catch DamnErr
        cd(curDir)
        error('Error: %s\n',DamnErr.message);
    end
end

% Prepare run of sesion subsampling

% Basic session information
Nsessions          = numel(GroupData.sessions);

% Samples are taken on a per subject basis, this grouping is used to
% ameliorate the temporal proximity sampling constraint i.e. I need to be
% able to subsample wavelet vectors in time, so I need to be able to space
% them out, but I also want a lot of samples. Therefore, I'm subsampling
% across 4 subjects, at 2 sessions each, instead of subsampling at the
% level of each session.

SesionsPerGroup  = 2;
L                = Nsessions/SesionsPerGroup;

SessionGroups    = cell(Nsessions/SesionsPerGroup,1);
SessionGroups{1} = [1 2];   % achilles
SessionGroups{2} = [3 4];   % buddy
SessionGroups{3} = [5 6];   % cicero
SessionGroups{4} = [7 8];   % gatsby

% Sampling parameters (following Berman "runEmbeddingSubSampling.m"
trainingSetSize    = parameters.trainingSetSize;  % for creating final embedding
training_numPoints = parameters.training_numPoints;
stepSz             = parameters.stepSz;

numPerDataSet      = round(trainingSetSize/L);
numModes           = GroupData.params.PCA.lastEigen;
numPeriods         = GroupData.params.wavelets.nbins;

trainingSetData    = zeros(numPerDataSet*L,numModes*numPeriods);

% set parameters for converting distance matrix to probability matrix
u            = parameters.training_perplexity;
tol          = GroupData.params.tSNE.tol;

% update the number of points per group to be used training to be the size
% of the smallest group's data after skipping by stepSz
GroupSubSampledSizes = nan(L,1);
for ii = 1:L
    tmp = 0;
    for jj = 1:SesionsPerGroup
        tmp = tmp + size(GroupData.sessions(SessionGroups{ii}(jj)).Data,1);
    end
    GroupSubSampledSizes(ii) = floor(tmp/stepSz);
end
training_numPoints            = min([GroupSubSampledSizes; training_numPoints]);
parameters.training_numPoints = training_numPoints;
fprintf(1,'Samples per Group of %d sessions: %d\n',SesionsPerGroup,training_numPoints);

%%
CUDAtSNE_dir = sprintf('%s%sMotion_Tracking%sBehaviorClassification%sUtils%scudaTSNE',...
    getenv('Dropbox'),filesep,filesep,filesep,filesep);

T1 = tic;
% for each session...
for kk = 1:L
    
    % Pooled data for this across group
    Amps        = cell(numel(SessionGroups{kk}),1);
    
    currentIdx = (1:numPerDataSet) + (kk-1)*numPerDataSet;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% START: Encapsulate: "file_embeddingSubSampling.m" %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for each session...
    for jj = 1:numel(SessionGroups{kk})
        
        % Pooled wavelet data for session
        Amps{jj}  = cat(1,GroupData.sessions(SessionGroups{kk}(jj)).wavelets.power);
        
    end
    
    % for each session grouping compose samples 
    data = cat(1,Amps{:});
    data = data(:,1:numModes*numPeriods);
    
    % IMPORTANT: unlike Berman, I use linspace to generate the desired
    % number of samples. This makes sure the samples are chosen
    % approximately uniformly across time, as the expense of the
    % inter-sample interval varying +/- 1 sample. Berman fixes the sampling
    % interval and as a consequence throws out samples from the beginning
    % of the session. This is fine for him because his subjects are not
    % engaged in a task, whereas mine are, and I don't want to bias my data
    % toward one end of the session.
    
    % keep the middle section of the samples
    dataResid  = floor((size(data,1)-training_numPoints*stepSz)/2);
    signalIdx  = dataResid+1:stepSz:dataResid+training_numPoints*stepSz;
    
    % normalize rows to 1
    amps        = sum(data,2);
    signalData = bsxfun(@rdivide,data(signalIdx,:),amps(signalIdx));
    signalAmps = amps(signalIdx);
    
    clear data amps
    
    fprintf(1,'\tConverting Distance into Probability Matrix...');
    T0 = tic;
    ProbMat = wavelets2TsneGPU(single(signalData'),single(u),single(tol));
    fprintf(1,'done in %5.2f seconds.\n',toc(T0));
    if any(isnan(ProbMat(:)))
        fprintf(1,'Error in wavelets2TSNEGPU\n');
        keyboard
    end
    
    % Compute tSNE
    cd(CUDAtSNE_dir)
    T0 = tic;
    yData = tsne_p(ProbMat, single(2), single(2));
    fprintf(1,'tSNE embedding computed in %5.2f seconds.\n',toc(T0));
    cd(curDir)
    
    figure, plot(yData(:,1),yData(:,2),'.k'), axis square
    
    if DEBUG_SAMPLES
        fprintf(1,'Having a look...\n');
        keyboard
    end
    
    clear ProbMat
    
    [trainingSetData(currentIdx,:),~] = ...
            findTemplatesFromData(signalData,yData,signalAmps,...
                                numPerDataSet,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% END: Encapsulate: "file_embeddingSubSampling.m" %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

fprintf(1,'Saving Training Set to %s...',computeFolder);
    save([computeFolder filesep 'tSNE_trainingSet.mat'],...
            'trainingSetData')
fprintf(1,'Created Training Set in %9.6f seconds.\n',toc(T1));

fprintf(1,'Saving parameters for run to %s...',computeFolder);
save([computeFolder filesep 'parameters_test.mat'],...
        'parameters')
fprintf(1,'done.\n');

%% 4) run TSNE
fprintf(1,'Running ComputeTSNE.m...\n');

%%% DEVELOPMENT PIPELINE %%%
try
    keep parameters trainingSetData
catch
    clear
end
curDir = pwd;

% directory for caching variables
computeFolder = [getenv('Kinect2') filesep 'Data\Data10.10.2014' filesep 'DataProcessing'];

if ~exist('parameters','var')
    fprintf(1,'Attempting to load parameters...\n');
    
    cd(computeFolder)
    try
        TestFile = 'parameters_test.mat';
        tmp      = dir(TestFile);
        fprintf(1,'Group Data parameters created at %s.\n',tmp.date);
        fprintf(1,'Importing parameters features...');
        load([computeFolder filesep TestFile],'parameters')
        fprintf(1,'done.\n\n');
        cd(curDir)
    catch DamnErr
        cd(curDir)
        error('Error: %s\n',DamnErr.message);
    end
end

if ~exist('trainingSetData','var')
    fprintf(1,'Attempting to load training data...\n');
    
    cd(computeFolder)
    try
        TestFile = 'tSNE_trainingSet.mat';
        tmp      = dir(TestFile);
        fprintf(1,'tSNEmap_trainingData created at %s.\n',tmp.date);
        fprintf(1,'Importing tSNEmap_trainingSetData...');
        load([computeFolder filesep TestFile],'trainingSetData')
        fprintf(1,'done.\n\n');
        cd(curDir)
    catch DamnErr
        cd(curDir)
        error('Error: %s\n',DamnErr.message);
    end
end
if size(trainingSetData,2) < size(trainingSetData,1)
    % tranpose to make row major (C/C++)
    trainingSetData = trainingSetData';
end

%% Set up relevant parameters
numModes   = parameters.pcaModes;
numPeriods = parameters.numPeriods;

Wsub = numModes*numPeriods;    % the number of features

fprintf(1,'Number of features used in clustering = %d\n',Wsub);

% Convert wavelet power into probability matrix
u   = parameters.perplexity;
tol = parameters.sigmaTolerance;

fprintf(1,'Converting Distance into Probability Matrix...');
T0 = tic;
ProbMat = wavelets2TsneGPU(single(trainingSetData),single(u),single(tol));
if any(isnan(ProbMat(:)))
    fprintf(1,'Error in wavelets2TSNEGPU\n');
    keyboard
end
fprintf(1,'done in %5.2f seconds.\n',toc(T0));


%% run TSNE
fprintf(1,'Running CUDA tSNE...\n');

% due to memory constraints, I'll just keep what's necessary
keep parameters ProbMat u tol computeFolder

% TODO: figure out why this crashed unless I am in the dir with the
% compiled files!
curDir       = pwd; 
CUDAtSNE_dir = sprintf('%s%sMotion_Tracking%sBehaviorClassification%sUtils%scudaTSNE',...
    getenv('Dropbox'),filesep,filesep,filesep,filesep);
cd(CUDAtSNE_dir)

T0 = tic;
yData = tsne_p(ProbMat, single(2), single(2));
fprintf(1,'tSNE embedding computed in %5.2f seconds.\n',toc(T0));
cd(curDir)

% flag for saving output tSNEmap
savetSNEmap = true;

if savetSNEmap
    
    % composing embedding data structure to be saved
    tSNEmap.data   = yData;
    params.u       = u;
    params.tol     = tol;
    tSNEmap.params = params;

    fprintf(1,'Saving tSNE map to %s...',computeFolder);
    save([computeFolder filesep 'tSNEmap_test.mat'],...
            'tSNEmap')
    fprintf(1,'done.\n');
end

%% estimate density using Berman code
%set up range
kdNeighbors = parameters.kdNeighbors;

maxY = ceil(max(abs(yData(:)))) + 1;

% set density smoothing parameter to median of k-nearest neighbors
% distance
NS = createns(yData);
[~,D] = knnsearch(NS,yData,'K',kdNeighbors+1);

sigma = median(D(:,kdNeighbors+1));

[xx,density] = findPointDensity(yData,sigma,501,[-maxY maxY]);
% Convert to pdf
density = density./sum(density(:));
density(density<1e-9) = 0;

% view method (for class implementation)
tSNEmap_view = true;
if tSNEmap_view
    % generate figure and view from above
    figure, surf(xx,xx,double(density),'linestyle','none'),view(0,90),
        axis([-maxY maxY -maxY maxY])
        axis square
        set(gca,'xtick',[-fix(maxY) 0 fix(maxY)],...
                'ytick',[-fix(maxY) 0 fix(maxY)])
    xlabel('z_1','FontName','Arial','FontSize',18)
    ylabel('z_2','FontName','Arial','FontSize',18)


    newmap = jet();       %change as desired, e.g., flag(256)

    minz = double(min(density(:)));
    if minz > 0
      disp('Your range is all above 0, no change');
    else
      maxz = double(max(density(:)));
      if maxz < 0
        disp('Your range is all below 0, no change');
      else
        ncol = size(newmap, 1);
        zratio = (0 - minz) ./ (maxz - minz);
        zpos = max( round(zratio * ncol), 1);  %closest non-zero
        newmap(zpos,:) = [1 1 1];   %set there to white
        colormap(newmap);           %activate it
      end
    end
end
