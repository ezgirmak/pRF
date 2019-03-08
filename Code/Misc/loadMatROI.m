function vertices = loadMatRoi(roiPath, labelindex)
    for i = 1:length(roiPath)
        load(roiPath{i});
        labels = find(ismember(wang2015labels, labelindex));
        indices =   ismember(wang2015,labels);
        vertices = find(indices);
    end
    


end