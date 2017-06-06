function img = removeborder(img, frac)
[r, c] = size(img);
mask = true(r,c);
mask(1:round(r*frac),:) = false;
mask(round(r*(1-frac)):end,:) = false;
mask(:,1:round(c*frac)) = false;
mask(:,round(c*(1-frac)):end) = false;
img(~mask)=NaN;
end
