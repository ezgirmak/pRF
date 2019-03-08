clear all close all


homeDr = '/SCRATCH/VISCOG_FOLDERS/';
dataDr = '/SCRATCH/VISCOG_FOLDERS/118225/pRF';

cd(dataDr)
load('118225_pRF_Retinotopy_Gaussian2D_2mm_01Dec2018_wang_angle_eccentricity.mat');
orig = collated;
load('118225_pRF_Retinotopy_Gaussian2D_2mm_filtered_01Dec2018_wang_angle_eccentricity.mat');
filt = collated;
i_orig = find([orig.pRF(:).didFit]); disp(['fitted orig = ', num2str(length(i_orig))]);

%% cc
f1 = figure(1);  clf; %title('orig diagnostics');
subplot(2,2,1); hist([orig.pRF(i_orig).corr]); title('corr');
subplot(2,2,2); hist([orig.pRF(i_orig).angle]);  title('angle');
subplot(2,2,3); hist([orig.pRF(i_orig).radius]);  title('ecc');
subplot(2,2,4); hist([orig.pRF(i_orig).sigma]); title('sigma')
saveas(f1, 'diagnostics_orig_wang.png');
f2 = figure(2); title('orig polar');clf
polarplot([orig.pRF(i_orig).angle], ...
    [orig.pRF(i_orig).radius], 'ko'); hold on
saveas(f2, 'polar_orig_wang.png');

return

%% fitting quality

i_filt = find([filt.pRF(:).didFit]); disp(['fitted filtered = ',num2str(length(i_filt))]);
i_both = intersect(i_orig, i_filt); disp(['fitted intersect = ',num2str(length(i_both))]);
f3 = figure(3); clf
subplot(2,2,1)
plot([orig.pRF(i_both).corr], [filt.pRF(i_both).corr], '.', 'MarkerSize', 5); title('fitting quality'); hold on
plot([0 1], [0 1], 'k--');
xlabel('original');
ylabel('filtered');

%% eccentricity biases
subplot(2,2,2)
plot([orig.pRF(i_both).radius], [filt.pRF(i_both).radius], '.', 'MarkerSize', 5);title('eccentricity bias'); hold on
plot([0 15], [0 15], 'k--');
xlabel('original');
ylabel('filtered');

%% angles biases
subplot(2,2,3)
plot([orig.pRF(i_both).angle], [filt.pRF(i_both).angle], '.', 'MarkerSize', 5); title('angle bias');hold on
plot([-8 8], [-8 8], 'k--');
xlabel('original');
ylabel('filtered');

%% sigma biases
subplot(2,2,4)
plot([orig.pRF(i_both).sigma], [filt.pRF(i_both).sigma], '.', 'MarkerSize', 5);title('sigma bias'); hold on
plot([0 4], [0 4], 'k--');
xlabel('original');
ylabel('filtered');
saveas(f3, 'comparison_wang.png');

%% distortions
f4 = figure(4); clf


for v =  1:length(i_both)
polarplot([orig.pRF(i_both(v)).angle filt.pRF(i_both(v)).angle], ...
    [orig.pRF(i_both(v)).radius filt.pRF(i_both(v)).radius], 'k');
hold on
polarplot([orig.pRF(i_both(v)).angle ], ...
    [orig.pRF(i_both(v)).radius] , 'o'); hold on

polarplot([filt.pRF(i_both(v)).angle ], ...
    [filt.pRF(i_both(v)).radius] , '>'); 
end
saveas(f4, 'comparison_polar_wang.png');

