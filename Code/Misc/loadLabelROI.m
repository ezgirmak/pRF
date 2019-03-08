function vertices = loadLabelROI(roiPath)

vertices = [];
for i = 1:length(roiPath)
    vertices = [vertices readLabels(roiPath{i})']; % read .label file, zero-based indexing
end

vertices = vertices +1;

end