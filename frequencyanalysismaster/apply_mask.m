function [masked_I] = apply_mask(M, I)
% This function takes a binary 2D mask (M)
% and applies it onto every 3rd dimension of I 
% The size of I in 1rst and 2nd dimesion need to be identical to the size
% of M


%% Check whether inputs are correct 
if size(M,1) == size(I,1) && size(M,2) == size(I,2) && islogical(M)
else 
    sprintf('Your mask does not fit your recording')
end

%% Create a mask of NaNs
mask = double(M);
mask(~M) = nan; 

%% Apply mask to the matrix
masked_I = I.*mask; 
end 