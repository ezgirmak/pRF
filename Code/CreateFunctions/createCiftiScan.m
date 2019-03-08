function [scan] = createCiftiScan(boldPath, roiPath, scanOpt)
% [scan] = createFreeSurferScan(boldPath, roiPath)
% Helpful function to createScan.m. Extracts the scan information from
% FreeSurfer .nii, .nii.gz, .mgz files (boldPath) within a given .label
% ROI. 
%
% Inputs: 
%   boldPath            Path to BOLD information (.nii, .nii.gz, .mgz), 
%                       string
%   roiPath             Path to ROI information (.label), if empty string, 
%                       will extract full brain time courses, string.
% Output:
%   scan                A structure containing the provided scan 
%                       infromation as fields:
%       boldFile        Name of the BOLD data file (.nii, .nii.gz, .mgz), 
%                       string
%       boldSize        Size of the BOLD data in [nVolumes x y z] format, 
%                       numeric
%       nVols           Number of volumes, numeric
%       TR              Scan TR, seconds
%       dur             Total scan duration, seconds
%       t               Time vector of the scan in TRs, seconds
%       voxID           Voxel index number, numeric
%       vtc             Voxel time course in [nVolumes nVox] format,
%                       numeric
%
% Note:
% - Dependencies: <a href="matlab: web('https://github.com/vistalab/vistasoft/tree/master/external/freesurfer')">mrVista/FreeSurfer</a>

% Written by Kelly Chang - July 19, 2017


%% Extract Scan Information from .nii, .nii.gz, .mgz (and .label)
% load data
wbcmd = 'wb_command';

bold = double(getfield(ciftiopen(boldPath,wbcmd),'cdata'));


if ~all(cellfun(@isempty, roiPath))
    if isfield(scanOpt, 'roifiletype') && strcmp(scanOpt.roifiletype, 'cifti');
        vertices = loadCiftiROI(roiPath, scanOpt.labelindex);
    elseif isfield(scanOpt, 'roifiletype') && strcmp(scanOpt.roifiletype, 'nifti');
          vertices = loadNiftiROI(roiPath);
    elseif isfield(scanOpt, 'roifiletype') && strcmp(scanOpt.roifiletype, 'mat');
           vertices = loadMatROI(roiPath,  scanOpt.labelindex);
    else
         vertices = loadLabelROI(roiPath);
    end
else
    vertices = 0:(size(tc,2)-1); % select all vertices
end

[~,file,ext] = fileparts(boldPath);
scan.TR = 1;
scan.nVols = size(bold, 1); % number of volumes in the scan

scan.dur = scan.nVols*scan.TR; % scan duration, seconds
scan.t = 0:scan.TR:(scan.dur-scan.TR); % time vector, seconds
scan.boldFile = [file ext]; % name of bold data file
scan.voxID = vertices; % vertices, zero-based indexing
scan.vtc = bold(vertices+1, :)'; % extract time course of vertices, +1 to zero-based 


