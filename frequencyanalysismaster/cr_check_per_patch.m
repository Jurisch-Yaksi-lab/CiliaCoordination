function [] = cr_check_per_patch(lmatrix, phase_patch,xfield, yfield, phase_pd, mean_phases, cv_phases, pos, PowerSpec, theta_all,magn_all, wave_length, CBF)
% This function plots the underlying variables that make up the wave
% direction and wavelength

%% Set some variables
binEdges = 0:(2*pi)/12:2*pi; % For the polar plot

%% Create a new folder where all png's will be saved
mkdir([CBF.targetP, 'patchresults\']);
cd([CBF.targetP, 'patchresults\']);

%% Loop over all patches in the recording

for iPatch = 1:CBF.NumPatches
    
    %% Prepare for the quiver plot
    [X,Y] = meshgrid(1:size(xfield(:,:,iPatch),2),1:size(xfield(:,:,iPatch),1)); % yes, arguments are swapped because of matlab: https://stackoverflow.com/questions/28418095
    
    % Normalize the gradient: We want unit vectors pointing at the gradient
    % direction.
    norm = magnitude(xfield(:,:,iPatch), yfield(:,:,iPatch));
    [xfield_norm, yfield_norm] = deal(xfield(:,:,iPatch)./norm, yfield(:,:,iPatch)./norm);
    
    %% Prepare for the wavelength histogram 
    lambdas = 2*pi./magn_all{iPatch}; 
    
    %% Plot all information in two figures
    set(figure, 'units','normalized','outerposition',[0 0 1 1])
    
    % Plot the entire segmentation
    ax1 = subplot(4,3,1);
    imagesc(lmatrix, 'AlphaData',~((lmatrix==0)));  colormap(ax1, 'jet')
    title('Segmentation'),set(gca,'YDir','normal');
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    xlabel('\mum'), ylabel('\mum'),axis image,
    
    % Display the patch that is examined here
    ax2 = subplot(4,3,2);
    imagesc((CBF.lmatrix == iPatch)*iPatch, 'AlphaData',~((CBF.lmatrix == iPatch)==0)); colormap(ax2, 'jet'),
    caxis(ax2,[0 max(CBF.lmatrix(:))]), title(sprintf('patch %d', iPatch))
    colormap jet,axis image
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    xlabel('\mum'), ylabel('\mum')
    
    % Display the heatmap of the contribution of the patch frequency
    ax3 = subplot(4,3,3);
    map = flipud(squeeze(PowerSpec(:,:,pos(iPatch))));
    borders = edge(lmatrix==iPatch, 'canny');
    imagesc(map, 'AlphaData', ~borders), colormap(ax3, 'hot'); c = colorbar; hold on,
    set(gca,'YDir','normal');
    title('heatmap');
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    xlabel('\mum'), ylabel('\mum');  axis image
    c.Label.String = 'spectral power';
    caxis(ax3,[0 1800]); map = [];
        
    % Display the phase for this patch
    ax4 = subplot(4,3,4);
    imagesc(phase_patch(:,:,iPatch), 'AlphaData', ~isnan(phase_patch(:,:,iPatch)));
    colormap(ax4, 'hsv'), caxis(ax4,[0, 2*pi])
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    set(gca,'YDir','normal');
    xlabel('\mum'), ylabel('\mum'); axis image
    title('phase map'); 
    
    % Display the arrows for this patch
    ax5 = subplot(4,3,5);
    quiver(X,Y,xfield_norm,yfield_norm,'color','k')
    xlabel('\mum'), ylabel('\mum')
    title('arrows'), axis image
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    
    % Plot the color wheel
    subplot(4,3,6); 
    create_color_wheelV2()
      
    % Display the direction per pixel
    ax7 = subplot(4,3,7);
    imagesc(theta_all{iPatch}, 'AlphaData', ~isnan(theta_all{iPatch})), colormap(ax7, 'hsv'), caxis(ax7,[0, 2*pi])
    xlabel('\mum'), ylabel('\mum')
    title('direction per pixel'), axis image
    set(gca,'YDir','normal');
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    
    % Display the directions in a polar plot
    subplot(4,3,8);
    ax8 = polarhistogram(theta_all{iPatch},binEdges);
    title(sprintf('mean_phase = %0.2f, sv = %0.2f',wrapTo2Pi(mean_phases(iPatch)), cv_phases(iPatch)));
    ax8.Parent.ThetaAxisUnits = 'radians';
    
    % Display the direction summary per patch
    ax9 = subplot(4,3,9);
    imagesc(phase_pd, 'AlphaData', (lmatrix == iPatch)), colormap(ax9, 'hsv'), caxis(ax9,[0, 2*pi])
    xlabel('\mum'), ylabel('\mum')
    title('direction per patch'), axis image
    set(gca,'YDir','normal');
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
  
    % Display the magnitude per pixel
    ax10 = subplot(4,3,10);
    imagesc(magn_all{iPatch}, 'AlphaData', (lmatrix == iPatch)),colormap(ax10, 'parula'),
    xlabel('\mum'), ylabel('\mum'),c = colorbar;  caxis(ax10,[0 1.8])
    c.Label.String = 'Magnitude in [rad/pixel]';
    title('Image gradient magnitude per patch');
    set(gca,'YDir','normal');axis image
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    
    % Display magnitude histogram
    subplot(4,3,11);
    % Plot a histogram of wavelengths.
    h = histfit(lambdas(lambdas<CBF.max_wavelength),CBF.bin_size, 'kernel'); hold on, h(1).FaceColor = [ 0.5843 0.8157 0.9882];
    % Identify the peak of the histogram and extract the value at the peak.
    [~,p] = max(h(2).YData); max_lambdas = h(2).XData(p);  
    % Draw a red line for the most common wavelength in the patch
    Y = ylim; plot([max_lambdas,max_lambdas], [0 Y(2)],'r'); 
    max_lambdas = []; % Clear this variable to avoid overlap
    title(sprintf('patch %d: lambda = %0.2f pixels & %0.2f µm', iPatch, max_lambdas, max_lambdas*CBF.spatres));
    xlabel('Wavelength in [pixels]'); ylabel('pixel count')
    
    % Display the wavelength summary
    ax12 = subplot(4,3,12);
    imagesc(wave_length, 'AlphaData', (lmatrix == iPatch)), colormap(ax12, 'jet'), caxis(ax12,[0 5])
    xlabel('\mum'), ylabel('\mum'), c=colorbar;
    title(sprintf('lambda = %0.2f µm', mean(wave_length(lmatrix == iPatch)))), axis image
    set(gca,'YDir','normal');
    set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
    set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
    c.Label.String = 'wavelength in [µm]';
    
    % Save the plot
    suptitle(sprintf('patch %d', iPatch))
    saveas(gcf, sprintf('results_patch%d.png',iPatch));

close all,
    
end


