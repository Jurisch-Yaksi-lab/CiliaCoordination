function [Mask] = create_mask(mask,signal_size)
% Create mask using image-analysis tools 

% % Convert the mask into a logical one
mask = ~isnan(mask);

% Find pixels that are connected to one another
B = bwconncomp(mask,8);

% Select only parts that are sufficiently large
MASK = cellfun('length',B.PixelIdxList)> signal_size; % See if 1000 works well
    
% only keep those components which are larger then minsize pixels
B.PixelIdxList(~MASK) = [];

% Prepare the new mask 
Mask = zeros(B.ImageSize); 
for i = 1:length(B.PixelIdxList)
    Mask(B.PixelIdxList{i}) = 1; 
end 

