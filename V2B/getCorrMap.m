function  getCorrMap(pRFPath, tempPath, savename)

        wbcmd = 'wb_command';
        %load prf
        load(pRFPath);
        %load roi
        cii = ciftiopen(tempPath,wbcmd);
        ciiData = double(cii.cdata);
        tmp = NaN(size(ciiData));        tmp = tmp(:);
        corr = cat(1,collated.pRF(:).corr);
        index = cat(1,collated.pRF(:).id);
        
        tmp(index) = corr;

        newcii = cii; newcii.cdata = tmp; ciftisave(newcii,savename,'wb_command');
   

       
        
  




end