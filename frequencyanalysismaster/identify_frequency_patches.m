function [lmatrix,complist] = identify_frequency_patches(I,CBF)
% This function returns frequency patches given an input image I
% that shows frequencies as determined for a cilia carpet. 
% This function assumes that the noise region is masked by NaN values.
% In addition, we require frequency patches to have minimum size of
% minsize.
%
% output
% lmatrix: matrix with entries indicating the different patches with labels
% complist: structure containing indices of the patches
% minsize: minimum size of a template used to identify connected patches
%
% 2018 Written by Stephan Bialonski

%% Preprocessing
% Determine set of different frequencies which are present in the data.
freq_set = unique(I(:));

complist.PixelIdxList = {}; 

fprintf('%s: Identify frequency patches\n',mfilename)

tic
for i = 1:length(freq_set)
    % determine connected components
    P = I==freq_set(i);
    c = bwconncomp(P,8);
    
    % check whether connected components are larger then minsize pixels
    mask = cellfun('length',c.PixelIdxList)>CBF.minsize;
    c.PixelIdxList(~mask) = [];
    
    % store these connected components
    count = length(complist.PixelIdxList);
    complist.PixelIdxList(count+1:count+length(c.PixelIdxList)) = c.PixelIdxList;
end
toc 

% complement information in complist for later processing
complist.Connectivity = c.Connectivity;
complist.ImageSize = c.ImageSize;
complist.NumObjects = length(complist.PixelIdxList);

% label patches
lmatrix = labelmatrix(complist);

% Fill the holes in those patches 
[lmatrix, complist] = bmf_fill_holes_patches(lmatrix, complist);


% save results
save(fullfile(CBF.targetP,[CBF.name,'_result_frequency_patches.mat']),'complist','-v7.3');

%% visualize results
fprintf('plotting and saving...\n\n')

figure, imagesc(lmatrix, 'AlphaData', ~lmatrix == 0);
colormap jet, set(gcf,'color','w')
title('frequency patches as obtained by segmentation')
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
xlabel('\mum'), ylabel('\mum'); axis image

saveas(gcf, [CBF.targetP, CBF.name, '_figure_determine_frequency_patches.png']);
print( '-painters', [CBF.targetP, CBF.name, '_figure_determine_frequency_patches'], '-depsc');