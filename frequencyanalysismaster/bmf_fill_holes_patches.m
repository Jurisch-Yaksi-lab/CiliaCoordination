function [lmatrix_new, complist] = bmf_fill_holes_patches(lmatrix, complist)

% Count the number of patches
nPatch = unique(lmatrix); nPatch = nPatch(2:end);

% Prepare a new matrix
lmatrix_new = zeros(size(lmatrix));

for i_patch = 1:length(nPatch)
    
    % binary mask for this patch
    mask = (lmatrix==i_patch);
    
    % Resegment the patches
    % generate structuring element
    se = strel('diamond',10); % the higher the number, the more fine features will become blurred
    % dilate the mask image
    mask_dilated = imdilate(mask,se);
    % erode the mask again; thus size will be approximately preserved, but fine features and small gaps will be gone
    mask_eroded = imerode(mask_dilated,se);
    
    % Save the new patch segmentation
    lmatrix_new(mask_eroded) = i_patch;
end

% This loop needs to be separate since some patches are now overlapping
for i_patch = 1:length(nPatch)
    
    % Save the new indices
    complist.PixelIdxList{i_patch} = find(lmatrix_new == i_patch);
    
end

end