% demopRF.m

clear all; close all;

%% Variables

% subject information
subject = '156334';

spath = '/mnt/viscog/FineLab/V2B/156334/'
% directories
paths = createPaths(); % initialize paths structure
paths.data = spath; % path to demostration data directory
paths.results = fullfile(spath, 'pRF'); % path to output results directory
paths = createPaths(paths); % create paths if they do not already exist

% model options
opt = createOpt('Retinotopy');
opt.roi = 'ROIs.1.60.nii.gz';
opt.estHRF = 3; 
opt.parallel = true;

% %% Convert DemoVTC.mat to DemoData.vtc
% % This section exists because of GitHub's file size limit. This section can
% % be deleted or commented 'DemoData.vtc' has been created. 
% 
% if ~exist(fullfile(paths.data, 'DemoData.vtc'), 'file')
%     load(fullfile(paths.data, 'DemoVTC.mat')); % load .mat version of .vtc data
%     tmp = BVQXfile('new:vtc'); % ~!~BVQXTools Dependency~!~
%     flds = fieldnames(vtc);
%     for i = 1:length(flds)
%         tmp.(flds{i}) = vtc.(flds{i});
%     end
%     tmp.saveAs(fullfile(paths.data, 'DemoData.vtc'));
% end

%% Variables (cont.)

% scan options
scanOpt.boldPath = fullfile(paths.data, '/functional/MNINonLinear/Results/tfMRI_RETCW_7T_PA/tfMRI_RETCW_7T_PA.nii.gz'); %
scanOpt.matPath = fullfile(paths.data, '/apertures/apertures/RETCWsmall.mat');
scanOpt.roiPath = fullfile(paths.data, '/structural/MNINonLinear/fsaverage_LR59k/156334.L.BA.59k_fs_LR.label.gii');
% scanOpt.paradigm.x = 'stimImg';

scan = createScan(scanOpt, opt); % creating 'scan' structure

% seed options
seedOpt.mu = linspace(2, 4, 21);
seedOpt.sigma = linspace(1, 4, 20);
seedOpt.exp = 0.5; % fixed parameter

seeds = createSeeds(seedOpt); % creating 'seeds' structure

% hrf options
hrfOpt.type = 'BoyntonHRF';
hrfOpt.dt = scan(1).TR;

hrf = createHRF(hrfOpt); % creating 'hrf' structure

%% Estimate and Save pRFs

[collated] = estpRF(scan, seeds, hrf, opt);

saveName = fullfile(paths.results, ...
    sprintf('%s_pRF_%s_%s_%s', upper1(subject), opt.map, opt.model, ...
    datestr(now, 'ddmmmyyyy')));

safeSave(saveName, 'collated');

%% Visualizations

% estimated paramaters
paramOpt.params = {'mu', 'sigma', 'exp', 'corr'};
paramOpt.measure = @median;
plotParams(collated, paramOpt);

% predicted time course
plotPredicted(collated);

% fitted hrf
plotHRF(collated.hrf);
