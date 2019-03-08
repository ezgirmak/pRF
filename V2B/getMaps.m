function  getMaps(type, collated, tempPath, savename)
wbcmd = 'wb_command';
%load prf
% load(pRFPath);
%load roi
cii = ciftiopen(tempPath,wbcmd);
ciiData = double(cii.cdata);
tmp = NaN(size(ciiData));        tmp = tmp(:);
switch type
    case 'corr'
        
        
        corr = cat(1,collated.pRF(:).corr);
        index = cat(1,collated.pRF(:).id);
        
        tmp(index) = corr;
        
        
    case 'angle'
        
        ang = cat(1,collated.pRF(:).angle);
        index = cat(1,collated.pRF(:).id);
        
        tmp(index) = ang;
        
    case 'radius'
        
        rad = cat(1,collated.pRF(:).radius);
        index = cat(1,collated.pRF(:).id);
        
        tmp(index) = rad;
        
        
    case 'sigma'
        
        
        sigma = cat(1,collated.pRF(:).corr);
        index = cat(1,collated.pRF(:).id);
        
        tmp(index) = sigma;
        
end

savename = [savename '_'  type  '.dscalar.nii'];
newcii = cii; newcii.cdata = tmp; ciftisave(newcii,savename,'wb_command');


end