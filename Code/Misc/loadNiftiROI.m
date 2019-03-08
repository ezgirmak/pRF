function vertices = loadNiftiROI(roiPath, labelindex)

 for i = 1:length(roiPath)
    tmp = MRIread(roiPath{1}); 
    tmp = tmp.vol(:);
%   vertices(ismember(tmp,labelindex)) = 1;
    indices =   ismember(tmp,labelindex);
    vertices = find(indices);
 end

end