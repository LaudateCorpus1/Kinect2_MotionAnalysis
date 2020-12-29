function [kmData, PCs, lastEigen] = Kinect2PCs(kmData, Inspect)
%Kinect2PCs : 
% This code takes in a kmData structure containing a matrix representing 
% the time-series of variables [Ntimepoints x Nvariables], and partitioned 
% into epochs. It then performs PCA jointly upon all epochs to denoise the 
% signal. The output is a structure analogous to the inputs structure 
% with the data now being PC weights of the user defined compoments of 
% interest in descending order by eigen value magnitude.
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% kmData (required): a multi-structure, one per file, with fields:
%   .Data: a [NFrames x Nnodes(xyz)] matrix of kinematic data
%   .NFrames: the number of frames in the kinect2 data file
%   .fileInfo: a string detailing which file was loaded
%
% Inspect: a flag indicating whether or not to output a diagnostic figure
%   detailing the choice of the value for the second output argument.
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% kmData: a multi-structure, one per file, with fields:
%   .Data: A [NFrames x numel(Vars)] PC weights associated with frames
%   .NFrames: the number of frames in the kinect2 data file
%   .fileInfo: a string detailing which file was loaded
%
% PCs: a structure with the fields:
%   .meanData: the means of the original data
%   .stdData:  the standard deviations of the original data
%   .projM:    the matrix of PCs ordered by descending eignvalues
%   .eigV:     a vector of eigenvalues sorted in decreasing order
%
% lastEigen: a scaler index indicating the last significant eigen value
%            relative to the shuffle control

%%%%%%%%%%%%%%%%%%%%
%%% INPUT CHECKS %%%
%%%%%%%%%%%%%%%%%%%%
mfname = mfilename;
if nargin < 1
    error('%s: Input Error: No inputs supplied: printing help file...\n%s\n',mfname, help(mfname));
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

% Flag for whether or not to output PC figure
if nargin < 2
    Inspect = false;
end

%%% Transform the data using PCA %%%
NEpochs = numel(kmData);
kmPC    = cat(1,kmData(:).Data);

% normalize data
% 1) zero mean and unit std
PCs.meanData = mean(kmPC);
PCs.stdData  = std(kmPC);
kmPC         = kmPC - repmat(PCs.meanData,[size(kmPC,1) 1]);
kmPC         = kmPC./repmat(PCs.stdData,[size(kmPC,1) 1]);

% 2) all variabes [-0.5 0.5] with zero mean
% PCs.meanData = nan;
% PCs.stdData  = nan;
% kmPC = kmPC - min(kmPC(:));
% kmPC = kmPC/max(kmPC(:));
% kmPC = bsxfun(@minus,kmPC,mean(kmPC,1));

if size(kmPC, 2) < size(kmPC, 1)
    C = kmPC' * kmPC;
else
    C = (1 / size(kmPC, 1)) * (kmPC * kmPC');
end

[projM, lambda] = eig(C);
[lambda, ind] = sort(diag(lambda), 'descend');
projM           = projM(:,ind);
if ~(size(kmPC, 2) < size(kmPC, 1))
    projM = bsxfun(@times, kmPC' * projM, (1 ./ sqrt(size(kmPC, 1) .* lambda))');
end

% Compare eigenvalue series from data to that of shuffled surrogates
tmp = kmPC;
for ii = 1:size(kmPC,2)
    inds      = randperm(size(kmPC,1),size(kmPC,1));
    tmp(:,ii) = tmp(inds,ii);
end

C            = tmp' * tmp;
[~, lambda0] = eig(C);
lambda0      = sort(diag(lambda0), 'descend');

% output the index of the last PC greater than shuffle noise floor
lastEigen = find(lambda > lambda0,1,'last');

if Inspect    
    ExplainedVar = 100*(sum(lambda(1:lastEigen))./sum(lambda));
    
    figure('units','pixels','position',[400 150 800 600]),
    semilogy(lambda0,'b','linewidth',2)
    hold on, semilogy(lambda,'r','linewidth',2),
    semilogy(lastEigen, lambda(lastEigen),'og','markerfacecolor','g')
    
    text(0.625*length(lambda), 0.625*lambda(1),sprintf('Explained Variance: %4.1f',ExplainedVar));
    axis tight
    xlabel('Ordered PCs')
    ylabel('Eigen Value (a.u.)')
    title('Sorted Eigen Values Series: Data(red) vs. Shuffled Surrogates(blue)')
    
    set(gca,'FontSize',11)
    
end

% project onto principle components
kmPC       = kmPC * projM;

% Store eigenvectors and eigenvalues
PCs.projM  = projM;
PCs.eigV   = lambda;

% deal back into trials, accounting for nans...
counter = 1;

for ii = 1:NEpochs
    kmData(ii).Data = kmPC(counter:counter + kmData(ii).NFrames-1,:);
    counter         = counter + kmData(ii).NFrames;
end

end