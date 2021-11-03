function [max_wavelength, bin_size] = cr_patchwise_analysisV2(phase_patch, pos, PowerSpec, CBF)
% In this function, a gradient map is calculated for a static phase. 
% From this gradient map, both the wave direction and wavelength will be
% calculated and plotted. 
% Credits go to Stephan Bialonski for writing the codes to determine the
% wavelenght and wave direction. I have just adapted it to calculate them
% for a single patch at the time. 

% phase_patch = the phase map calculated at the peak frequency of that
% patch
% pos = the frequencies per patch in power spectral position 
% PowerSpec = the powerspectra 
% CBF = A structure containing necessary variables. 


%% Variables 
max_wavelength = 200; 
bin_size = 30; 
NumPatches = size(phase_patch,3); 

% Define this anonymous function 
magnitude = @(x,y) sqrt(x.^2 + y.^2);

close all, 

%% Calculate the wave direction per patch

[xfield, yfield, phase_pd, phase_sd] = wave_direction(CBF.lmatrix); %#ok<ASGLU>

% Save figures
set(figure(1), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(1), [CBF.targetP, CBF.name, '_figure_polarplot_patch_direction.png']);
set(figure(2), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(2), [CBF.targetP, CBF.name, '_figure_per_patch_direction.png']);

% Vector map 
vis_gradient()

%% Calculate the wavelength per patch

[wave_length] = wavelength(CBF.lmatrix);

% Save figures
set(figure(8), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(8), [CBF.targetP, CBF.name, '_figure_histogram_fit_patch_wavelength.png']);

save([CBF.targetP, CBF.name, '_results_wave_direction_and_lengths.mat'], 'wave_length', 'phase_pd', 'phase_sd');
close all, 

%% Nested functions
    function [xfield, yfield, phase_pd, phase_sd] = wave_direction(lmatrix)
        
        % Allocate for speed
        mean_phases = nan(1,NumPatches); % Circular mean per patch
        cv_phases  = nan(1,NumPatches); % Circular variance per patch
        phase_pd = nan(size(lmatrix)); % To store the new patch image that encodes mean directions
        phase_sd = nan(size(lmatrix)); % To store the new patch image that encodes the variance of directions
        binEdges = 0:(2*pi)/12:2*pi; % For the polar plot
        xfield = nan([size(lmatrix),NumPatches]);
        yfield = nan([size(lmatrix),NumPatches]);
        
       
        for i = 1:NumPatches
            
            % Determine directions of gradient vectors
            mask_patch = lmatrix == i;
            [theta_patch, xfield(:,:,i), yfield(:,:,i)] = determine_theta(phase_patch(:,:,i), 1);
            
            % Restrict directions to the frequency patch we want to investigate
            theta_patch = theta_patch(~isnan(theta_patch));
            
            % Determine mean phase angle
            if ~isempty(theta_patch)
                mean_phases(i) = angle(sum(exp(1i.*theta_patch(:))));
%                 max_phases(i) = angle(exp(1i.*theta_patch(:)));
                cv_phases(i) = 1-abs(mean(exp(1i.*theta_patch(:)))); % determine circular standard deviation of mean phase angle
            else
                mean_phases(i) = nan;
                cv_phases(i) = nan;
            end
            
            % Store the mean direction in a frequency patch image
            phase_pd(mask_patch) = mean_phases(i)+pi; % shift the phases to [0, 2pi] interval
            phase_sd(mask_patch) = cv_phases(i);
            
            % Plot polarplots and directions per patch to check for
            figure(1), subplot(7,6,i), polarhistogram(theta_patch+pi,binEdges), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));
            figure(2), subplot(7,6,i), imagesc(phase_pd, 'AlphaData', mask_patch), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum'); axis image
            colormap hsv,  caxis([0 2*pi])
        end
        
        %% Create figures
        
        % Phase direction map
        figure,
        imagesc(phase_pd, 'AlphaData', ~isnan(phase_pd))
        colormap hsv, hold on
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image,  caxis([0 2*pi])
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_direction.png']);       
%         print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_direction'], '-depsc');
        
        
        % Phase direction variance map
        figure,
        imagesc(phase_sd, 'AlphaData', ~isnan(phase_pd))
        colormap gray,  c = colorbar;
        hold on, caxis([0 1]), c.Label.String = 'standard deviation';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_direction_variance.png']); 
        
        % Phase per patch map
        figure,
        imagesc(nansum(phase_patch,3), 'AlphaData', ~isnan(phase_pd))
        colormap hsv,  c = colorbar; caxis([0 2*pi])
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        hold on
        c.Label.String = 'phase angle';
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_phasemap.png']); 
        print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_phasemap'], '-depsc');
        
        % Phase direction and variance map
        figure,
        imagesc(phase_pd, 'AlphaData', 1 - phase_sd)
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        colormap hsv,  c = colorbar;  caxis([0 2*pi])
        c.Label.String = 'wave direction'; axis image
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_combined_direction_and_variance.png']); 
        print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_combined_direction_and_variance'], '-depsc');
        
        % Heatmap per patch
        figure('units','normalized','outerposition',[0 0 1 1])
        for n = 1:NumPatches
            map(:,:,n) = squeeze(PowerSpec(:,:,pos(n)));
            
            % Plot the heatmap
            borders = edge(lmatrix==n, 'canny');
            subplot(7,6,n), imagesc(squeeze(map(:,:,n)), 'AlphaData', ~borders), colormap hot; c = colorbar; hold on,
            title(sprintf('heatmap for patch %i', n))
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum');  axis image
            c.Label.String = 'spectral power';
            caxis([0 1800])
        end
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_heatmap.png']);
        
       
    end

    function [wave_length] = wavelength(lmatrix)
        %% Determine gradient magnitudes for the frequency patch
        
        % Allocate memory
        max_lambdas = nan(1, NumPatches);
        magn_rad = nan(size(lmatrix)); 
        magn_pix = nan(size(lmatrix)); 
        wave_length = nan(size(lmatrix));
        
        for i = 1:NumPatches
            
            % Determine magnitudes
            magn = magnitude(squeeze(xfield(:,:,i)), squeeze(yfield(:,:,i))); % [rad/pixel]
            
            % Determine metachronal wave lengths
            magn_patch = magn(~isnan(magn));
            lambdas = 2*pi./magn_patch; % Wavelength in [pixels] for all pixels
            mask_patch = ~isnan(magn);
            
            % Store the mean direction in a frequency patch image
            magn_rad(mask_patch) = magn_patch; % [rad/pixel]
            magn_pix(mask_patch) = lambdas; % [pixel]
            
            % Attempt to identify one common wavelength per patch (sometimes a patch is empty)
            try
                
                % Plot the histogram of vectors for every patch
                figure(8), subplot(7,6,i),  
                h = histfit(lambdas(lambdas<max_wavelength),bin_size, 'kernel'); % Plot the histogram with kernel of wavelengths
                [~,p] = max(h(2).YData); max_lambdas(i) = h(2).XData(p); hold on, % Identify the peak of the histogram and extract the value at the peak.
                Y = ylim; plot([max_lambdas(i),max_lambdas(i)], [0 Y(2)]); % Draw a red line for the most common wavelength in the patch
                title(sprintf('patch %d: lambda %0.2f', i, max_lambdas(i))); xlabel('Wavelength in [pixels]'),
                wave_length(lmatrix==i) = max_lambdas(i)*CBF.spatres; % wavelength in [µm] per patch
                
            catch
                sprintf('No lambdas found in patch %d', i);
            end
                                
        end
        figure(8),subtitle('Histogram of magnitudes per patch');
        
        %% Create figures
        
        % Magnitude in [rad/pixels]
        figure, imagesc(magn_rad, 'AlphaData', ~isnan(magn_rad));
        c = colorbar;
        c.Label.String = 'Magnitude in [rad/pixel]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum')
        axis image
        hold on
        caxis([0 1.8])
        title('Image gradient magnitude per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_magnitude_rad_pixel.png']);
        
        % Magnitude in [pixels]
        figure, imagesc(magn_pix, 'AlphaData', ~isnan(magn_pix));
        c = colorbar;
        c.Label.String = 'Magnitude in [pixel]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 500])
        title('Image gradient magnitude per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_magnitude_pixel.png']);
        
        % Wavelength in micrometer per patch
        figure, imagesc(wave_length, 'AlphaData', ~(lmatrix==0));
        c = colorbar; colormap jet;
        c.Label.String = 'Wavelength [µm]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 10])
        title('wavelength per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_wavelength_per_patch.png']);
        print( '-painters', [CBF.targetP, CBF.name, '_figure_wavelength'], '-depsc');
               
    end

    function [] = vis_gradient()
        
        % Visualize the gradient field
        Xfield = nansum(xfield,3);
        Yfield = nansum(yfield,3);
        
        % Create 2D support points (i.e. a grid) where the velocity vectors are
        % located.
        [X,Y] = meshgrid(1:size(Xfield,2),1:size(Yfield,1)); % yes, arguments are swapped because of matlab: https://stackoverflow.com/questions/28418095
        
        % Normalize the gradient: We want unit vectors pointing at the gradient
        % direction.
        norm = magnitude(Xfield, Yfield);
        [xfield_norm, yfield_norm] = deal(Xfield./norm, Yfield./norm);
        
        figure,
        imagesc(nansum(phase_patch,3), 'AlphaData', ~(nansum(phase_patch,3)) == 0)
        colormap parula;  c = colorbar;
        hold on
        quiver(X,Y,xfield_norm,yfield_norm,'color','w')
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum')
        axis image  
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_gradient.png']);
        
    end

end