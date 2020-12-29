function k2data2d = loadKinect2dataRot(kinect2BinFile)
%loadKinect2data2d :
% This routine reads in Kinect2 kinematic model Quaternion data contained 
% within a binary file. This data can be used for extracting joint angles.
%
%%%%%%%%%%%%%%
%%% Inputs %%%
%%%%%%%%%%%%%%
% kinect2BinFile:   a binary file containing kinematic data formatted according
%                   to the Kinect2Acq specificaitons:
%                       25 joints, 2 coordinates, float
%
%%%%%%%%%%%%%%%
%%% Outputs %%%
%%%%%%%%%%%%%%%
% k2dataRot:     a structure array with fields:
%         .Data: a [Nnodes(wxyz) x NFrames] matrix of quaternion data
%       .Nnodes: the number of nodes
%      .Ncoords: the number of coordinates for each node
%      .Nframes: the number of frames in this file
%     .fileInfo: a string indicating the file input by the user
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com

mfname = mfilename;
if nargin < 1
    error('%s: Input Error: No inputs supplied: printing help file...\n%s\n',mfname, help([mfname '.m']));
end

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

if ~strcmp(filename(end-3:end),'_Rot')
    error('%s: The file ''%s'' is not a Quaternion binary file.', mfname, kinect2BinFile);
end

% Load Kinect2 data
fprintf(1,'Loading %s...',filename);
fid         = fopen(kinect2BinFile,'r');
Nnodes      = 25;
Ncoords     = 4;
SizeOfFloat = 4;
NFrames     = fileInfo.bytes/(Nnodes*Ncoords*SizeOfFloat);
fprintf(1,'done\n');

K2Data      = fread(fid,[Ncoords*Nnodes NFrames],'float');
fclose(fid);

% check correctness of file
if (NFrames - floor(NFrames)) > 0
    error(['%s: Incorrect File: The computed number of frames in ''%s'' is '...
           'not an integer value.'], mfname, kinect2BinFile);
end

% Compose output
k2data2d.Data     = K2Data;
k2data2d.Nnodes   = Nnodes;
k2data2d.Ncoords  = Ncoords;
k2data2d.Nframes  = NFrames;
k2data2d.fileInfo = kinect2BinFile;
