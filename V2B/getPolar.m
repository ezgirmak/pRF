% Code to get eccentricity and angle from cartesian coordinates.
function collated = getPolar(data)
%load data
%load(filepath);
tmp = data;
x  = [tmp.pRF(:).xMu];
y = [tmp.pRF(:).yMu];

[angle radius] = cart2pol(x,y);
angle = num2cell(angle);
radius = num2cell(radius);
tmp.pRF = rmfield(tmp.pRF, {'xMu', 'yMu'});
[tmp.pRF.angle] = angle{:};
[tmp.pRF.radius] = radius{:};

collated = tmp;
end


