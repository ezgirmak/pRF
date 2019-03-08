function [scan] = createScan(scanOpt, opt)
% [scan] = createScan(scanOpt, opt)
%
% Creates a structure 'scan' containing information about the scan(s) given
% by the corresponding 'scanOpt.boldPath' and 'scanOpt.matPath' through
% one of two methods
%
% METHOD 1: PARADIGM
% Will create a 'scan' structure based on the given
% 'scanOpt.paradigm.<funcOf>' sequence(s), will create a stimulus image
% from the given paradigm sequence(s)
%
% METHOD 2: STIMULUS IMAGE
% Will create a 'scan' structre based on pre-defined stimulus image(s)
%
% Inputs:
%   scanOpt                  A structure containing option to create the
%                            'scan' structure with fields:
%       boldPath             Path(s) to all BrainVoyager or FreeSurfer BOLD
%                            files(s), string
%       matPath              Path(s) to all .mat paradigm files, string
%       roiPath              Path(s) to all BrainVoyager or FreeSurfer ROI
%                            file(s), string
% -------------------------------------------------------------------------
% METHOD 1: PARADIGM
%       paradigm             A structure containing variable name(s) of the
%                            paradigm sequence(s) located within the
%                            paradigm files; optional, specifying will
%                            let the code create the stimulus image
%            <funcOf         Stimulus paradigm sequence of the field
%             parameters>    specified 'funcOf' parameter name
%
% METHOD 2: STIMULUS IMAGE
%       stimImg              Name of the stimulus image variable within the
%                            paradigm file(s), string; optional, specifying
%                            will use this variable as the stimulus image
%       dt                   Time step size, numeric (default:
%                            scan.dur/size(stimImg,1))
%       order                Order of the stimulus image dimensions
%                            (default: [nVolumes <opt.model's funcOf>])
%       funcOf               A structure containing the stimulus function
%                            of dimension range as fields
% -------------------------------------------------------------------------
%   opt                      A structure containing option for pRF model
%                            fitting containing fields:
%       model                Model name, also the function name to be
%                            fitted string
%       roi                  Name(s) of the ROI files if fitting within
%                            ROI(s), string
%
% Output:
%   scan                     A structure with length 1xN where N is the
%                            length of 'scanOpt.boldPath' containing the
%                            .vtc and scan's information with fields:
%       matFile              Name of .mat file , string
%                            (i.e., 'Subj1_Paradigm_Set1.mat')
%       paradigm             A structure containing the paradigm sequences
%                            for each stimulus dimension given as fields:
%           <funcOf          Stimulus paradigm sequence, should be given in
%             parameters>    units that are to be estimated for each
%                            stimulus dimension, blanks should be coded as
%                            NaNs, numeric
%       k                    A structure containing the unique stimulus
%                            values for each stimulus dimension given as
%                            fields:
%           <funcOf          Unique stimulus values for each stimulus
%             parameters>    dimension, excludes NaNs
%       funcOf               A structure containing the actual function of
%                            parameters as matrices scan given the model as
%                            fields:
%           <funcOf          Full function of stimulus values for each
%             parameters>    stimulus dimension, meshgrid applied if
%                            multiple funcOf parameters
%       boldFile             Name of the BOLD file, string
%                            (i.e., 'Subj1_Set1.<vtc/nii>')
%       boldSize             Size of the functional data
%       nVols                Number of volumes in the scan
%       dur                  Total scan duration, seconds
%       TR                   TR of the scan, seconds
%       dt                   Time step for the paradigm or stimulus image,
%                            seconds
%       t                    Time vector of the scan in TRs, seconds
%       voxID                Voxel index number
%       vtc                  Voxel time course
%       stimImg              A M x Ni x ... x Nn matrix where M is the
%                            number of volumes of the scan and Ni through
%                            Nn is the length(scan.paradigm.<funcOf>) or
%                            the desired resolution of the stimulus image
%                            for each stimulus dimension
%
% Note:
% - Dependencies: <a href="matlab:
% web('http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools/232-getting-started.html')">BVQXTools/NeuroElf</a>, <a href="matlab: web('https://github.com/vistalab/vistasoft/tree/master/external/freesurfer')">mrVista/FreeSurfer</a>

% Written by Kelly Chang - June 23, 2016
% Edited by Kelly Chang - September 1, 2017
%Edited by Ezgi Yucel, Mar 7, 2019

%% Input Control

if ~isfield(scanOpt, 'boldPath') || isempty(scanOpt.boldPath)
    error('No bold files selected');
end

if ischar(scanOpt.boldPath)
    scanOpt.boldPath = {scanOpt.boldPath};
end

if ~isfield(scanOpt, 'matPath') && isempty(scanOpt.matPath)
    error('No .mat files selected');
end

if ischar(scanOpt.matPath)
    scanOpt.matPath = {scanOpt.matPath};
end

if length(scanOpt.boldPath) ~= length(scanOpt.matPath)
    error('All bold files must have corresponding .mat files');
end

if ~isfield(scanOpt, 'roiPath')
    scanOpt.roiPath = {''};
end

if ischar(scanOpt.roiPath)
    scanOpt.roiPath = {scanOpt.roiPath};
end

if isempty(opt.roi) && ~all(cellfun(@isempty, scanOpt.roiPath))
    error('No ''opt.roi'' when ''scanOpt.roiPath'' is specified');
end

if ~isempty(opt.roi) && all(cellfun(@isempty, scanOpt.roiPath))
    error('No ''scanOpt.roiPath'' when ''opt.roi'' is specified');
end

if isfield(scanOpt, 'paradigm') && isfield(scanOpt, 'stimImg')
    error('Cannot specify both ''scanOpt.paradigm.<var>'' and ''scanOpt.stimImg''');
end

if isfield(scanOpt, 'paradigm') && ~isstruct(scanOpt.paradigm)
    error('Must specify variable name(s) for ''scanOpt.paradigm.<var>''');
end

paramNames = eval(opt.model);
if isfield(scanOpt, 'paradigm') && ~all(ismember(paramNames.funcOf, fieldnames(scanOpt.paradigm)))
    errFlds = setdiff(paramNames.funcOf, fieldnames(scanOpt.paradigm));
    error('Must specify paradigm.<var> for variable(s): %s', strjoin(errFlds, ', '));
end

if isfield(scanOpt, 'paradigm') && any(structfun(@isempty, scanOpt.paradigm))
    errFlds = fieldnames(scanOpt.paradigm);
    error('Must specify ''paradigm.<var>'' variable name(s) for: %s', ...
        strjoin(errFlds(structfun(@isempty, scanOpt.paradigm)), ', '));
end

if isfield(scanOpt, 'stimImg') && isempty(scanOpt.stimImg)
    error('Must specify a variable name for ''scanOpt.stimImg''');
end

if isfield(scanOpt, 'stimImg') && ~isfield(scanOpt, 'order')
    scanOpt.order = ['nVols' paramNames.funcOf];
end

%% Bold File Name(s)

[~,file,ext] = cellfun(@fileparts, scanOpt.boldPath, 'UniformOutput', false);
boldFile = strcat(file, ext);

%% Creating 'scan' Structure

flds = fieldnames(scanOpt);
stimImgMethod = char(flds(ismember(flds, {'paradigm', 'stimImg'})));
[~,~,software] = cellfun(@fileparts, scanOpt.boldPath, 'UniformOutput', false);
for i = 1:length(scanOpt.boldPath)
    if ~opt.quiet
        fprintf('Loading: %s\n', boldFile{i});
    end
    
    tmp= switchSoftware(scanOpt.boldPath{i}, scanOpt.roiPath, scanOpt);
   
    %     scan = tmp; % save output
    switch stimImgMethod % load stimulus data
        case 'paradigm' % specifying with paradigm sequence
            scan(i) = createStimImg(tmp, scanOpt, i, opt);
        case 'stimImg' % extracting pre-made stimImg
            scan(i) = extractStimImg(tmp, scanOpt, i, opt);
    end
end
end


function scan = switchSoftware(boldPath, roiPath, scanOpt)
% scan = switchSoftware(boldPath, roiPath, scanOpt)
% 
% Helpful local function to createScan.m. Extracts the scan information from
% given BOLD files for given ROI's.
%
% Inputs: 
%   boldPath            Path to BOLD information (.vtc), string
%   roiPath             Path to ROI information (.voi), if empty string, 
%                       will extract full brain time courses, string.
%   scanOpt             Scan options object
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
% Note:
% - Dependencies: <a href="matlab:
% web('http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools/232-getting-started.html')">BVQXTools/NeuroElf</a>
% <a href="matlab: web('https://github.com/vistalab/vistasoft/tree/master/external/freesurfer')">mrVista/FreeSurfer</a>
% <a href="matlab: web('https://www.humanconnectome.org/software/get-connectome-workbench.html')">Connectome Workbench/ HCP</a>
% <a href="matlab: web('https://www.artefact.tk/software/matlab/gifti/')">GiftiToolbox/ HCP</a>

% Written by Kelly Chang - July 19, 2017
% Edited by Ezgi Yucel, Mar 7, 2019


switch scanOpt.boldfiletype % loading bold data
    case 'BV' % BrainVoyager
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
        
        
        scan.boldSize = size(bold.VTCData); % size of the .vtc data
        scan.nVols = bold.NrOfVolumes; % number of volumes in the scan
        scan.TR = bold.TR/1000; % seconds
        scan.voxID = [vtc.id]; % voxel id number (linearized)
        scan.vtc = [vtc.vtcData]; % voxel time course
    case 'HCP'% FreeSurfer or HCP
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
        
        scan.TR = 1;
        scan.nVols = size(bold, 1); % number of volumes in the scan
        scan.voxID = vertices; % vertices, zero-based indexing
        scan.vtc = bold(vertices+1, :)'; % extract time course of vertices, +1 to zero-based
        
    case 'FS' %freesurfer
        bold = MRIread(boldPath); % load .nii file
        tc = squeeze(bold.vol);
        tc = permute(tc, [ndims(tc) 1:(ndims(tc)-1)]); % bold time course
        scan.boldSize = size(tc); % size of the bold data
        tc = reshape(tc, size(tc,1), []);
        
        if ~all(cellfun(@isempty, roiPath))
            if isfield(scanOpt, 'roifiletype') && strcmp(scanOpt.roifiletype, 'cifti')
                vertices = loadCiftiROI(roiPath);
            elseif isfield(scanOpt, 'roifiletype') && strcmp(scanOpt.roifiletype, 'nifti')
                vertices = loadNiftiROI(scanOpt.roiPath, scanOpt.labelindex);
            else
                vertices = loadLabelROI(roiPath);
            end
        else
            vertices = 1:(size(tc,2)); % select all vertices
        end
        
        scan.nVols = bold.nframes; % number of volumes in the scan
        scan.TR = bold.tr/1000; % seconds
        scan.voxID = vertices; % vertices, zero-based indexing
        scan.vtc = tc(:,vertices);
    otherwise
        error('Unrecognized filetype: %s', scanOpt.boldfiletype);
end
[~,file,ext] = fileparts(boldPath);
scan.boldFile = [file ext]; % name of bold data file
scan.dur = scan.nVols*scan.TR; % scan duration, seconds
scan.t = 0:scan.TR:(scan.dur-scan.TR); % time vector, seconds
end