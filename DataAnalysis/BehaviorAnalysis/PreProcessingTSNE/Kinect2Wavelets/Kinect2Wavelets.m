function waveletFeatures = Kinect2Wavelets(kmData, varargin)
%Kinect2Wavelets :
%
% This code takes in a matrix representing time-series of kinematic variables,
% and outputs a matrix of wavelet coefficients for each kinematic variable 
% i.e. It is a map from kinematic variables to wavelet spectra. This is a one-to-many
% mapping. This routine is written with Kinect/Kinect2 data in mind.
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% kmData (required): a multi-structure, one per file, with fields:
%   .Data: a [NFrames x Nvars] matrix of kinematic data
%       -this data may be pre-processed e.g. PCA has been run
%   .NFrames: the number of frames in the kinect2 data file
%   .fileInfo: a string detailing which file was loaded
%
% varagin: string "key/value" pairs for overriding default parameters
%
%%% Free Parameters %%%
% Wavelet decomposition Free Parameters (see defaults below)
% Fs:       Samples per second
% fMin:     Minimum frequency for wavelet
% fMax:     Maximum frequency for wavelet
% nbins:    number of requested wavelets
% omega:    non-dimensional parameter of Morlet wavelet
% Inspect:  flag for producing a diagonostic figure
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% waveletFeatures: a structure with fields:
%   .wavelets: a multi-structure, one per epoch, with fields:
%       .power: an [Nsamples x Nwavelets*Nvars] matrix e.g. for 25 wavelets and
%          5, 3D nodes, there would be 375 columns, 75 for each node (xyz).
%       .ts:  the timestamps, in seconds for each frame in an epoch
%       
%   .params: a structure with fields:
%       .(free parameters)
%       .freqs: the frequencies associated with the wavelets
%       .scales: morlet wavelet scaling factors
%       .fileInfo: a cell of strings detailing which files was loaded

%%%%%%%%%%%%%%%%%%%%
%%% Input checks %%%
%%%%%%%%%%%%%%%%%%%%
mfname = mfilename;
if nargin < 1
    error('%s: Input Error: No inputs supplied: printing help file...\n%s\n',mfname, help([mfname '.m']));
end

if ~isfield(kmData,{'Data'; 'NFrames'; 'fileInfo'})
    error('%s: InputError: field error: printing help file...\n%s\n',mfname, help(mfname));
end

% check if there are more rows than columns in each file. Transpose as
% necessary and warn
for ii = 1:numel(kmData)
    if size(kmData(ii).Data,1) < size(kmData(ii).Data,2)
        kmData(ii).Data = kmData(ii).Data';
        warning('%s: kmData(%d).Data transposed.',mfname,ii);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Free Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet decomposition Free Parameters
Defaults.Fs      = 30;
Defaults.fMin    = 0.25;
Defaults.fMax    = Defaults.Fs;
Defaults.nbins   = 25;
Defaults.omega   = 5;
Defaults.Inspect = false;

%%%%%%%%%%%%%%%%%%%%%%%
%%% PARSE ARGUMENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
    % There is no argument validation herein.
    Defaults = parseArgs(Defaults,varargin);
end
% return parameters stored in Defaults
v2struct(Defaults);

% log parameters
Params.Fs    = Fs;
Params.fMin  = fMin;
Params.fMax  = fMax;
Params.nbins = nbins;
Params.omega = omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Wavelet decomposition of kinematic model Frames %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output
Nepochs = numel(kmData);
waveletFeatures = struct('wavelets',cell(1,1),'params',cell(1,1));
    Wavelets    = struct(...
                         'power',   cell(Nepochs,1),...
                         'ts',      cell(Nepochs,1));

fprintf(1,'Calculating wavelets...\n');
T0 = tic;
reshapeTmp = cell(size(kmData(1).Data,2),1);

Files = cell(Nepochs,1);
for ii = 1:Nepochs
    fprintf(1,'\tEpoch %d out of %d\n',ii,Nepochs);
    % access data
    data      = kmData(ii).Data;
    Files{ii} = kmData(ii).fileInfo;
    
    % having a look at the power spectrum
    if ii == 1
        [ amp, ~, freqs, scales, t] = getWavelet(data,Fs, fMin, fMax, nbins, omega);
        Params.freqs  = freqs;
        Params.scales = scales;
    else
        [ amp, ~, ~, ~, t] = getWavelet(data,Fs, fMin, fMax, nbins, omega);
    end
    
    % reshape wavelet power to [ nSamples x nWavelets ]
    for jj = 1:size(amp,3)
        reshapeTmp{jj} = amp(:,:,jj)';
    end
    amp  = cat(2,reshapeTmp{:});
    % log output
    Wavelets(ii).power   = amp;
    Wavelets(ii).ts      = t;
end
Params.fileInfo = Files;

% compose top-level structure
waveletFeatures.wavelets = Wavelets;
waveletFeatures.params   = Params;
fprintf(1,'...done in %6.3f seconds.\n',toc(T0));

% Diagnostic figure: inspect the first 3 wavelet Features
if Inspect
    % the number of level sets for the contour plots
    NlevelSets = 20;
    
    for ii = 1:Nepochs
        
        ts       = waveletFeatures.wavelets(ii).ts;
        freqs    = waveletFeatures.params.freqs;
        features = length(freqs);
        
        tmp = waveletFeatures.wavelets(ii).power(:,1:3*features);
        
        % Assemble freqs and determine axis scaling
        w   = sqrt(tmp);
        w0  = w(:,1:features);
        w1  = w(:,features+1:2*features);
        w2  = w(:,2*features+1:3*features);
        
        %scale data to max across per feature dimensions
        MaxW = max([w0(:) w1(:) w2(:)]);
        w0   = w0./MaxW(1);
        w1   = w1./MaxW(2);
        w2   = w2./MaxW(3);
        
        tmpFig = figure('Name',sprintf('Epoch %d',ii),...
            'Position',get(0,'ScreenSize'),'color','w');
        
        ax1 = subplot(3,1,1);
        [~, H] =contourf(ts,freqs,w0',NlevelSets);axis xy
        set(H,'LineStyle','none');
        title('1st wavelet feature'),ylabel('Freq (s^{1/2})')
        axis tight
        set(gca,'FontName','Arial','FontSize',18)
        
        ax2 = subplot(3,1,2);
        [~, H] =contourf(ts,freqs,w1',NlevelSets);axis xy
        set(H,'LineStyle','none');
        title('2nd wavelet feature'),ylabel('Freq (s^{1/2})')
        axis tight
        set(gca,'FontName','Arial','FontSize',18)

        ax3 = subplot(3,1,3);
        [~, H] =contourf(ts,freqs,w2',NlevelSets);axis xy
        set(H,'LineStyle','none');
        title('3rd wavelet feature'),ylabel('Freq (s^{1/2})'),
        xlabel('Time (seconds)')
        axis tight
        set(gca,'FontName','Arial','FontSize',18)

        h = get(get(gcf,'children'),'xlabel');
        set(cat(1,h{:}),'FontSize',18)
        h = get(get(gcf,'children'),'ylabel');
        set(cat(1,h{:}),'FontSize',18)
        h = get(get(gcf,'children'),'title');
        set(cat(1,h{:}),'FontSize',18)

        linkaxes([ax1 ax2 ax3],'xy')
        pause
        close(tmpFig)
    end
end

end