clear all; close all


homeDr = '/mnt/viscog/FineLab/V2B/';
dataDr = '/mnt/viscog/FineLab/V2B/412528/pRF';

load('412528_pRF_Retinotopy_Gaussian2D_2mm_29Nov2018_angle_eccentricity.mat');
orig = collated;
i_orig = find([orig.pRF(:).didFit]); disp(['fitted orig = ', num2str(length(i_orig))]);


cd(dataDr)
load('412528_pRF_Retinotopy_Gaussian2D_2mm_filtered_30Nov2018_angle_eccentricity.mat');
filt = collated;
i_filt = find([filt.pRF(:).didFit]); disp(['fitted filt = ', num2str(length(i_orig))]);

i_both = intersect(i_orig, i_filt); disp(['fitted intersect = ',num2str(length(i_both))]);




%% plot data
alist = linspace(-pi, pi, 180);
rlist = linspace(0, 11.4, 100);


for a=1:length(alist)-1
    for r=1:length(rlist)-1
    orig_ind = find([orig.pRF.angle]>alist(a) & [orig.pRF.angle]<=alist(a)+1 ...
        & [orig.pRF.radius]>alist(a) & [orig.pRF.angle]<=alist(a)+1);
    cval_orig(a,r)=mean([orig.pRF(orig_ind).corr]);
     filt_ind = find([filt.pRF.angle]>alist(a) & [filt.pRF.angle]<=alist(a)+1 ...
        & [filt.pRF.radius]>alist(a) & [filt.pRF.angle]<=alist(a)+1);
    cval_filt(a,r)=mean([filt.pRF(filt_ind).corr]);
    
    end
end

figure(1); clf
subplot(1,2,1)
image(cval_orig.*256)
subplot(1,2,2)
image(cval_orig.*256)


