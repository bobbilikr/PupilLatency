function m = meanDisk(img, xc, yc, r)

[sy sx] = size(img);
xmin = max(1, floor(xc-r));
xmax = min(sx, ceil(xc+r));
ymin = max(1, floor(yc-r));
ymax = min(sy, ceil(yc+r));
img = img(ymin:ymax, xmin:xmax); 
xc = xc - xmin + 1;
yc = yc - ymin + 1;


[x y] = meshgrid(1:size(img,2), 1:size(img,1));
mask = (x-xc).^2 + (y-yc).^2 < r.^2;

% Compute mean
m = sum(sum(double(img) .* mask)) / sum(mask(:));

end

