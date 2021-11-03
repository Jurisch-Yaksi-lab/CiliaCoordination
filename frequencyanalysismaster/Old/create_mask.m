function [mask] = create_mask(mask)
% Create mask using image-analysis tools 

% Convert the mask into a logical one
mask = ~isnan(mask);

% basic morphological operations
se = strel('diamond',1); % remove isolated pixels
mask=imerode(mask,se);

se = strel('disk',1); % dilate/erode/dilate
for k=1:2
    mask=imdilate(mask,se);
end
mask=imfill(mask,'holes');
for k=1:4
    mask=imerode(mask,se);
end
for k=1:2
    mask=imdilate(mask,se);
end
