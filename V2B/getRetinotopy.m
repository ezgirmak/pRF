function  getRetinotopy(scanOpt)
 for i = 1:length(scanOpt.roiPath)
         load(scanOpt.roiPath{i});
        %load prf
        load(pRFPath);
        %load roi
        roi_nii = MRIread(roiPath);

        corr = cat(1,collated.pRF(:).corr);
        index = cat(1,collated.pRF(:).id);
        corr_nii = roi_nii;

        tmp = NaN(size(roi_nii.vol));

        tmp = tmp(:);
        tmp(index) = corr;
        corr_nii.vol = reshape(tmp, [size(roi_nii.vol)]);

        saveName = fullfile(paths.data, '156334_pRF_Retinotopy_Gaussian2D_19Nov2018_corr_map.nii.gz');
        MRIwrite(corr_nii,   saveName);
    end




end