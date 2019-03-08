
clear all; close all;

%add related paths
%% Variables
% subject information
subject = '146432';
wbcmd = 'wb_command';

% directories
root = getenv('ROOT');
addpath(genpath(root));
addpath(genpath([root '/HCPpipelines-3.27.0/global/matlab/gifti-1.6/']))
paths = createPaths(); % initialize paths structure
paths.main = pwd();
paths.data = fullfile(root, subject, '7T_RET_2mm_preproc/MNINonLinear/Results'); % path to data directory
paths.results = fullfile(root, subject, 'pRF'); % path to output results directory

% model options
opt = createOpt('Retinotopy');
opt.model = 'Gaussian2D';
opt.roi ='';% fullfile(paths.data, '3T_Structural_preproc_extended/T1w/156334/label' );
opt.estHRF = 0; % use default HRF for now, change it to 2 when you want perfection
opt.parallel = true;
opt.roi = true;
opt.corrThr = 0.25;
opt.freeList = {'-8<xMu<8', '-8<yMu<8', '0.1<sigma<3'};
opt.cmd = wbcmd; %unix Caret7command used in various functions
%% Variables (cont.)

% scan options


scanOpt.boldfiletype = 'HCP'; %BV if brainvoyager, HCP if HCP, FS is freesurfer
scanOpt.datatype = '2mm';
scanOpt.roiPath = fullfile(root, 'atlas.mat' );
scanOpt.roifiletype = 'mat';
scanOpt.labelindex ={'V1v','V1d','V2v','V2d','V3v','V3d', 'LO1','hV4', 'V3B','V3A'}; %Wang Atlas
%{'V1','V2','V3','V4'};  % for cifti files, where the number refers to (for example) the Brodmann area, for nifti files number refers to freesurfer color lookuptable
% %https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT

scanOpt.stimImg = 'stim';
scanOpt.dt = 1;

% scanStimNames = {'RETBARsmall.mat',  'RETCCWsmall.mat',  'RETCONsmall.mat',  'RETCWsmall.mat',  'RETEXPsmall.mat', 'RETBARsmall.mat'};
% scanBoldDir =   {'tfMRI_RETBAR1_7T_AP', 'tfMRI_RETCCW_7T_AP', 'tfMRI_RETCON_7T_PA',  'tfMRI_RETCW_7T_PA', ...
%                   'tfMRI_RETEXP_7T_AP',  'tfMRI_RETBAR2_7T_PA'};
% scanBoldFilenames = { 'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll.dtseries.nii ', ...
%                      'tfMRI_RETCCW_7T_AP_Atlas_MSMAll.dtseries.nii',...
%                      'tfMRI_RETCON_7T_PA_Atlas_MSMAll.dtseries.nii', ...
%                      'tfMRI_RETCW_7T_PA_Atlas_MSMAll.dtseries.nii', ...
%                      'tfMRI_RETEXP_7T_AP_Atlas_MSMAll.dtseries.nii',...
%                      'tfMRI_RETBAR2_7T_PA_Atlas_MSMAll.dtseries.nii'}\

switch scanOpt.datatype
    case '2mm'
        scanStimNames{1} = {'RETBARsmall.mat'};scanBoldDir{1} ={'tfMRI_RETBAR1_7T_AP'};
        scanBoldFilenames{1} = {'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll.dtseries.nii'};
        
        
        scanStimNames{2} = {'RETCCWsmall.mat'};scanBoldDir{2} ={'tfMRI_RETCCW_7T_AP'};
        scanBoldFilenames{2} = { 'tfMRI_RETCCW_7T_AP_Atlas_MSMAll.dtseries.nii'};
        
        scanStimNames{3} = {'RETCONsmall.mat'};scanBoldDir{3} ={'tfMRI_RETCON_7T_PA'};
        scanBoldFilenames{3} = { 'tfMRI_RETCON_7T_PA_Atlas_MSMAll.dtseries.nii'};
        
        scanStimNames{4} = {'RETCWsmall.mat'};scanBoldDir{4} ={'tfMRI_RETCW_7T_PA'};
        scanBoldFilenames{4} =  {'tfMRI_RETCW_7T_PA_Atlas_MSMAll.dtseries.nii'};
        
        scanStimNames{5} = {'RETEXPsmall.mat'};scanBoldDir{5} ={'tfMRI_RETEXP_7T_AP'};
        scanBoldFilenames{5} =  {'tfMRI_RETEXP_7T_AP_Atlas_MSMAll.dtseries.nii'};
        
        scanStimNames{6} = {'RETBARsmall.mat'};scanBoldDir{6} ={'tfMRI_RETBAR2_7T_PA'};
        scanBoldFilenames{6} =  {'tfMRI_RETBAR2_7T_PA_Atlas_MSMAll.dtseries.nii'};
    case '2mm_filtered'
        scanStimNames{1} = {'RETBARsmall.mat'};scanBoldDir{1} ={'tfMRI_RETBAR1_7T_AP'};
        scanBoldFilenames{1} = {'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_filtered.dtseries.nii'};
        
        
        scanStimNames{2} = {'RETCCWsmall.mat'};scanBoldDir{2} ={'tfMRI_RETCCW_7T_AP'};
        scanBoldFilenames{2} = { 'tfMRI_RETCCW_7T_AP_Atlas_MSMAll_filtered.dtseries.nii'};
        
        scanStimNames{3} = {'RETCONsmall.mat'};scanBoldDir{3} ={'tfMRI_RETCON_7T_PA'};
        scanBoldFilenames{3} = { 'tfMRI_RETCON_7T_PA_Atlas_MSMAll_filtered.dtseries.nii'};
        
        scanStimNames{4} = {'RETCWsmall.mat'};scanBoldDir{4} ={'tfMRI_RETCW_7T_PA'};
        scanBoldFilenames{4} =  {'tfMRI_RETCW_7T_PA_Atlas_MSMAll_filtered.dtseries.nii'};
        
        scanStimNames{5} = {'RETEXPsmall.mat'};scanBoldDir{5} ={'tfMRI_RETEXP_7T_AP'};
        scanBoldFilenames{5} =  {'tfMRI_RETEXP_7T_AP_Atlas_MSMAll_filtered.dtseries.nii'};
        
        scanStimNames{6} = {'RETBARsmall.mat'};scanBoldDir{6} ={'tfMRI_RETBAR2_7T_PA'};
        scanBoldFilenames{6} =  {'tfMRI_RETBAR2_7T_PA_Atlas_MSMAll_filtered.dtseries.nii'};
        
end


%specify the directories for the stimulus and timecourses
for scanNum = 1:length(scanStimNames)
    scanOpt.matPath = fullfile(root, 'apertures', scanStimNames{scanNum}); %this particular file is permute(stim, [3 1 2]) and then
    scanOpt.boldPath =   fullfile(paths.data,scanBoldDir{scanNum}, scanBoldFilenames{scanNum});
    
    scan(scanNum) = createScan(scanOpt, opt); % creating 'scan' structure
    x = linspace(-8,8,size(scan(scanNum).stimImg,3));
    y = linspace(-8,8, size(scan(scanNum).stimImg,2));
    [scan(scanNum).funcOf.x,scan(scanNum).funcOf.y] = meshgrid(x,y);
end

% seed options
seedOpt.xMu = linspace(-8, 8, 10);
seedOpt.yMu =  linspace(-8, 8, 10);
seedOpt.sigma = linspace(1, 4, 10);
seedOpt.exp = 1; % fixed parameter

seeds = createSeeds(seedOpt); % creating 'seeds' structure

% hrf options
hrfOpt.funcName = 'BoyntonHRF';
hrfOpt.dt = scan(1).TR;
hrfOpt.type = 'vision';

hrf = createHRF(hrfOpt); % creating 'hrf' structure

%% Estimate and Save pRFs

[collated] = estpRF(scan, seeds, hrf, opt);
saveName = fullfile(paths.results, ...
    sprintf('%s_pRF_%s_%s_%s_%s', upper1(subject), opt.map, opt.model, scanOpt.datatype, ...
    datestr(now, 'ddmmmyyyy')));


safeSave(saveName, 'collated', '-v7.3');

%% Visualizations

% estimated paramaters
paramOpt.params = {'xMu', 'yMu' , 'sigma','corr '};
paramOpt.measure = @median;

%get eccentricity and polar angle
collated = getPolar(collated);
saveName = fullfile(paths.results, ...
    sprintf('%s_pRF_%s_%s_%s_%s_angle_eccentricity', upper1(subject), opt.map, opt.model, scanOpt.datatype, ...
    datestr(now, 'ddmmmyyyy')));


safeSave(saveName, 'collated', '-v7.3');
tempPath = fullfile(basepath,subject,'3T_Structural_preproc/MNINonLinear/fsaverage_LR32k',[ num2str(subject) '.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii']);
%get Correlation maps
name  = fullfile(paths.results, sprintf('%s_pRF_%s_%s_%s_%s', upper1(subject), opt.map, opt.model, scanOpt.datatype, ...
    datestr(now, 'ddmmmyyyy')));
getMaps('angle', collated, tempPath, name);
getMaps('sigma', collated, tempPath, name);
getMaps('radius', collated, tempPath, name);
getMaps('corr', collated, tempPath, name);

