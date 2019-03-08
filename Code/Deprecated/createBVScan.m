function [scan] = createBVScan(boldPath, roiPath)
% [scan] = createBVScan(boldPath, roiPath)
% 
% Helpful function to createScan.m. Extracts the scan information from
% BrainVoyager .vtc files (boldPath) within a given .voi ROI.
%
% Inputs: 
%   boldPath            Path to BOLD information (.vtc), string
%   roiPath             Path to ROI information (.voi), if empty string, 
%                       will extract full brain time courses, string.
%
% Output:
%   scan                A structure containing the provided scan 
%                       infromation as fields:
%       boldFile        Name of the BOLD data file (.vtc), string
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
% - Dependencies: <a href="matlab:
% web('http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools/232-getting-started.html')">BVQXTools/NeuroElf</a>
% <a href="matlab: web('https://github.com/vistalab/vistasoft/tree/master/external/freesurfer')">mrVista/FreeSurfer</a>
% <a href="matlab: web('https://www.humanconnectome.org/software/get-connectome-workbench.html')">Connectome Workbench/ HCP</a>
% <a href="matlab: web('https://www.artefact.tk/software/matlab/gifti/')">GiftiToolbox/ HCP</a>

% Written by Kelly Chang - July 19, 2017
% Edited by Ezgi Yucel, Mar 7, 2019

%% Extract Scan Information from .vtc (and .voi)

bold = BVQXfile(boldPath); % load .vtc file
if ~all(cellfun(@isempty, roiPath))
    vtc = []; % initialize vtc
    for i = 1:length(roiPath)
        vtc = [vtc VTCinVOI(bold, BVQXfile(roiPath{i}))];
    end
else
    vtc = fullVTC(bold);
end

[~,file,ext] = fileparts(boldPath);
scan.boldFile = [file ext]; % name of bold data file
scan.boldSize = size(bold.VTCData); % size of the .vtc data
scan.nVols = bold.NrOfVolumes; % number of volumes in the scan
scan.TR = bold.TR/1000; % seconds
scan.dur = scan.nVols*scan.TR; % scan duration, seconds
scan.t = 0:scan.TR:(scan.dur-scan.TR); % time vector, seconds
scan.voxID = [vtc.id]; % voxel id number (linearized)
scan.vtc = [vtc.vtcData]; % voxel time course