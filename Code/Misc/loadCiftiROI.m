function vertices = loadCiftiROI(roiPath, labelindex)
wbcmd = 'wb_command';
vertices = [];
for j = 1:length(roiPath)
    
    roi =  double(getfield(ciftiopen(roiPath{j},wbcmd),'cdata'));
    grot=fileparts(roiPath{j});
    if (size(grot,1)==0)
        grot='.';
    end
    for i = 1:length(labelindex)
        
        tmp = [grot '/' labelindex{i} '.dscalar.nii'];
        tic
        unix([wbcmd ' -cifti-label-to-roi ' roiPath{j} ' ' tmp ' -name ' labelindex{i} ]);
        toc
        logical = double(getfield(ciftiopen(tmp,wbcmd),'cdata'));
        indices = ismember(logical,1);
        vertices = cat(1,vertices, find(indices)) ;
    end
    
end

%     for i = 1:length(roiPath)
%         tmp =  double(getfield(ciftiopen(roiPath{i},wbcmd),'cdata'));
%         indices =   ismember(tmp,labelindex);
%         vertices = find(indices);
%     end
%


end