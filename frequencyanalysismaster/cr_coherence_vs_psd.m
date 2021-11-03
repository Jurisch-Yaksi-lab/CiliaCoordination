function [map] = cr_coherence_vs_psd(val,  PeakPos, PowerSpec, CBF)
% This function goes with cr_coherence_ref
% Take the previous reference pixels and scatterplot the coherence versus PSD for that given pixel  

% OBS: I need to save the ?? 


%% Heatmaps of reference pixels 
map = zeros(CBF.x, CBF.y, length(CBF.ref));

for iref = 1:size(CBF.ref,1)
    r = CBF.ref(iref,:);
    map(:,:,iref) = heatmap_coherence(); 
end

%   Save variables
save([CBF.targetP, CBF.name, '_result_val_map_scatter'], 'val', 'map');


%% Scatterplot

for iref = 1:size(CBF.ref,1)
    
    % Define the reference pixels
    r = CBF.ref(iref,:);
    xx = r(1);
    yy = r(2);
    
    x_val = val(:,:,iref); x_val = x_val(~isnan(CBF.mask));
    y_map = map(:,:,iref); y_map = y_map(~isnan(CBF.mask));
    
    figure,
    scatter(x_val,y_map,6, 'filled', 'k', 'MarkerFaceAlpha', 0.3);
    ylim([0 1800]);
    ylabel('Frequency contribution'); xlabel('Coherence')
    title(sprintf('Coherence versus frequency contribution: pixel [%d,%d]',xx,yy))
    
    %   Save the figure
    saveas(gcf, [CBF.targetP,'ref_pix_coherence\' sprintf('pixel_%d_%d_scatter.png', xx,yy)]);
    
    close all,
    
end

%% Nested function definitions

    function [map] = heatmap_coherence()
        
        % Define the reference pixels
        xx = r(1);
        yy = r(2);
        
        % Collect the relative contribution of one pixel's peak to all other pixels
        pos = PeakPos(yy,xx);
        map = squeeze(PowerSpec(:,:,pos));
        
        % Plot the heatmap
        figure, imagesc(map, 'AlphaData', CBF.mask), colormap hot; c = colorbar;
        freq = pos * (CBF.Fs / CBF.nframe); 
        title(sprintf('heatmap for [%d,%d], f = %0.1f', xx,yy, freq))
        c.Label.String = 'spectral power';
        caxis([0 1800])
        
        % Save the figure
        saveas(gcf, [CBF.targetP,'ref_pix_coherence\',sprintf('pixel_%d_%d_heatmap.png', xx,yy)]);
        close all,
        
    end


end