function [scan] = createHCPScan(boldPath)
% [scan] = createHCPScan(boldPath)
%10.28.18
% Now includes support for CIFTI files as well 
%
% Helpful function to createScan.m. Extracts the scan information from
% HCP CIFTI files (boldPath) 
%
% Inputs: 
%   boldPath            Path to BOLD information (.nii, .nii.gz, .mgz), 
%                       string
%
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
%                * workbench
%                * HCPpipelines-3.27.0   
% Edited by Ezgi Yucel - Nov 05, 2018
% Written by Kelly Chang - July 19, 2017


%% Extract Scan Information from .nii, .nii.gz, .mgz (and .label)
% load data
wbcmd = 'wb_command';

bold = double(getfield(ciftiopen(boldPath,wbcmd),'cdata'));



%kendrick kay data manipulations
if ~iscell(data)
  data = {data};
end

% calc
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numvxs = prod(xyzsize);

%bold = MRIread(boldPath); % load .nii file
tc = squeeze(bold.vol);
tc = permute(tc, [ndims(tc) 1:(ndims(tc)-1)]); % bold time course
scan.boldSize = size(tc); % size of the bold data
tc = reshape(tc, size(tc,1), []);

if ~all(cellfun(@isempty, roiPath))
    vertices = [];
    for i = 1:length(roiPath)
        vertices = [vertices readLabels(roiPath{i})']; % read .lable file, zero-based indexing
    end
else
    vertices = 0:(size(tc,2)-1); % select all vertices
end

[~,file,ext] = fileparts(boldPath);
scan.boldFile = [file ext]; % name of bold data file
scan.nVols = bold.nframes; % number of volumes in the scan
scan.TR = bold.tr/1000; % seconds
scan.dur = scan.nVols*scan.TR; % scan duration, seconds
scan.t = 0:scan.TR:(scan.dur-scan.TR); % time vector, seconds
scan.voxID = vertices; % vertices, zero-based indexing
scan.vtc = tc(:,vertices+1); % extract time course of vertices, +1 to zero-based indexing