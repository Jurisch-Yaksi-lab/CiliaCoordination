function [] = cr_patchwise_analysis(phase_patch, pos, PowerSpec, CBF)
% cr_patchwise_analysis(phase_patch, pos, Fs, nframe, targetP, lmatrix, PowerSpec)

% Use the frequency segments to calculate a phase for that patch
% Use this phasemap to calculate the gradient and subsequently the wave
% direction as well as wavelength. 


% In this function, a gradient map is calculated for a static phase. 
% From this gradient map, both the wave direction and wavelength will be
% calculated. 


%% Variables 
% minsize = 400; 
NumPatches = size(phase_patch,3); 

% Define this anonymous function 
magnitude = @(x,y) sqrt(x.^2 + y.^2);

close all, 

%% Calculate the wave direction per patch

[xfield, yfield] = wave_direction(CBF.lmatrix);

% Save figures
set(figure(1), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(1), [CBF.targetP, CBF.name, '_figure_polarplot_patch_direction.png']);
set(figure(2), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(2), [CBF.targetP, CBF.name, '_figure_per_patch_direction.png']);

% Vector map 
vis_gradient()
set(figure(13), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(13), [CBF.targetP, CBF.name, '_figure_gradient_per_patch.png']);

%% Calculate the wavelength per patch

wavelength(CBF.lmatrix);

% Save figures
set(figure(8), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(8), [CBF.targetP, CBF.name, '_figure_histogram_fit_patch_wavelength.png']);

set(figure(9), 'units','normalized','outerposition',[0 0 1 1])
saveas(figure(9), [CBF.targetP, CBF.name, '_figure_magnitude_per_patch.png']);

close all, 

%% Nested functions
    function [xfield, yfield] = wave_direction(lmatrix)
        
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
                cv_phases(i) = 1-abs(mean(exp(1i.*theta_patch(:)))); % determine circular standard deviation of mean phase angle
            else
                mean_phases(i) = nan;
                cv_phases(i) = nan;
            end
            
            % Store the mean direction in a frequency patch image
            phase_pd(mask_patch) = mean_phases(i)+pi; % shift the phases to [0, 2pi] interval
            phase_sd(mask_patch) = cv_phases(i);
            
            % Plot polarplots and directions per patch to check for
            figure(1), subplot(6,5,i), polarhistogram(theta_patch+pi,binEdges), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));
            figure(2), subplot(6,5,i), imagesc(phase_pd, 'AlphaData', mask_patch), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum'); axis image
            colormap hsv
        end
        
        %% Create figures
        
        % Phase direction map
        figure,
        imagesc(phase_pd, 'AlphaData', ~isnan(phase_pd))
        colormap hsv, hold on
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_direction.png']);       
        print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_direction'], '-depsc');
        
        
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
        colormap hsv,  c = colorbar;
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        hold on
        c.Label.String = 'phase angle';
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_phasemap.png']); 
        
        % Phase direction and variance map
        figure,
        imagesc(phase_pd, 'AlphaData', 1 - phase_sd)
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        colormap hsv,  c = colorbar;
        c.Label.String = 'wave direction'; axis image
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_combined_direction_and_variance.png']); 
        
        % Heatmap per patch
        figure('units','normalized','outerposition',[0 0 1 1])
        for n = 1:NumPatches
            map(:,:,n) = squeeze(PowerSpec(:,:,pos(n)));
            
            % Plot the heatmap
            borders = edge(lmatrix==n, 'canny');
            subplot(6,5,n), imagesc(squeeze(map(:,:,n)), 'AlphaData', ~borders), colormap hot; c = colorbar; hold on,
            title(sprintf('heatmap for patch %i', n))
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum');  axis image
            c.Label.String = 'spectral power';
            caxis([0 1800])
        end
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_heatmap.png']);
        
       
    end

    function [] = wavelength(lmatrix)
        %% Determine gradient magnitudes for the frequency patch
        
        % Allocate memory
        max_lambdas = nan(1, NumPatches);
        magn_rad = nan(size(lmatrix)); % Do something here
        magn_pix = nan(size(lmatrix)); % Do something here
        wave_fwhm = nan(size(lmatrix));
        wave_length = nan(size(lmatrix));
        
        for i = 1:NumPatches
            % Determine magnitudes
            magn = magnitude(squeeze(xfield(:,:,i)), squeeze(yfield(:,:,i))); % Should be [rad/pixel]?
            
            % determine metachronal wave lengths
            magn_patch = magn(~isnan(magn));
            lambdas = 2*pi./magn_patch; % this works! wavelength in [pixels]
            mask_patch = ~isnan(magn);
            
            % Store the mean direction in a frequency patch image
            magn_rad(mask_patch) = magn_patch; % [rad/pixel]
            magn_pix(mask_patch) = lambdas; % [pixel]
            
            try
                % Plot the histogram of vectors for every patch
                figure(8), subplot(6,5,i),  h = histfit(lambdas(lambdas<200),30, 'kernel'); % !!!IMPORTANT!!!
                [v,p] = max(h(2).YData);
                max_lambdas(i) = h(2).XData(p); hold on,
                
                % Compute Full-Width-At-Half-Maximum [= width of the histogram]
                critval=v/2;
                iw1=find(h(2).YData>critval,1,'first');
                iw2=find(h(2).YData>critval,1,'last');
                if ~isempty(iw1) && ~isempty(iw2)
                    FWHM(i)=(h(2).XData(iw2)-h(2).XData(iw1)) / length(lambdas(lambdas<200)); % [pixels]
                else
                    FWHM(i)=nan;
                end
                
                Y = ylim; plot([max_lambdas(i),max_lambdas(i)], [0 Y(2)]); % Draw a line for the selected wavelenght
                plot([h(2).XData(iw1),h(2).XData(iw2)],[critval, critval], 'kx-' ); % Draw a line for the FWHM.
                title(sprintf('patch %d: lambda %0.2f FWHM %0.2f', i, max_lambdas(i), FWHM(i)))
                xlabel('Wavelength in [pixels]'),
                wave_fwhm(lmatrix==i) = FWHM(i);
                wave_length(lmatrix==i) = max_lambdas(i);
            catch
                sprintf('No lambdas found in patch %d', i);
            end
            
            
            
            % Plot the magnitude per patch
            figure(9), subplot(6,5,i), imagesc(magn_rad, 'AlphaData', mask_patch);
            c = colorbar; caxis([0 2]), hold on
            c.Label.String = 'Magnitude in [rad]';
            title('Image gradient magnitude per patch');
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum')
            axis image
        end
        figure(8),subtitle('Histogram of magnitudes per patch');
        figure(9),subtitle('Magnitudes per patch');
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
        figure, imagesc(wave_length*CBF.spatres, 'AlphaData', ~(lmatrix==0));
        c = colorbar; colormap jet;
        c.Label.String = 'Wavelength [µm]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 10])
        title('wavelength per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_wavelength_per_patch.png']);
        print( '-painters', [CBF.targetP, CBF.name, '_figure_wavelength_per_patch'], '-depsc');
        
        % Normalized FWHM per patch
        figure, imagesc(wave_fwhm, 'AlphaData', ~(lmatrix==0));
        c = colorbar;
        c.Label.String = 'FWHM';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 2])
        title('FWHM of the wavelength histogram per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_fmwh.png']); 
        
        % Wavelength in micrometer per patch adjusted for FWHM
        figure, imagesc(wave_length*CBF.spatres, 'AlphaData', 1 - wave_fwhm);
        c = colorbar;
        c.Label.String = 'Magnitude in [µm]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 30])
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_combined_wavelength_and_fmwh.png']); 
        
    end

    function [] = vis_gradient()
        %% Visualize the gradient field
        
        % All in one plot
        %         Xfield = nansum(xfield,3);
        %         Yfield = nansum(yfield,3);
        %
        %         % Create 2D support points (i.e. a grid) where the velocity vectors are
        %         % located.
        %         [X,Y] = meshgrid(1:size(Xfield,2),1:size(Yfield,1)); % yes, arguments are swapped because of matlab: https://stackoverflow.com/questions/28418095
        %
        %         % Normalize the gradient: We want unit vectors pointing at the gradient
        %         % direction.
        %         norm = magnitude(Xfield, Yfield);
        %         [xfield_norm, yfield_norm] = deal(Xfield./norm, Yfield./norm);
        %
        %         figure,
        %         imagesc(nansum(phase_patch,3), 'AlphaData', ~(nansum(phase_patch,3)) == 0)
        %         colormap parula;  c = colorbar;
        %         hold on
        %         quiver(X,Y,xfield_norm,yfield_norm,'color','w')
        
        % Plot per patch
        for i = 1:NumPatches
            
            % Create 2D support points (i.e. a grid) where the velocity vectors are
            % located.
            [X,Y] = meshgrid(1:size(xfield(:,:,i),2),1:size(yfield(:,:,i),1)); % yes, arguments are swapped because of matlab: https://stackoverflow.com/questions/28418095
            
            % Normalize the gradient: We want unit vectors pointing at the gradient
            % direction.
            norm = magnitude(xfield(:,:,i), yfield(:,:,i));
            [xfield_norm, yfield_norm] = deal(xfield(:,:,i)./norm, yfield(:,:,i)./norm);
            
            figure(13), subplot(6,5,i)
            imagesc(squeeze(phase_patch(:,:,i)), 'AlphaData', ~isnan(squeeze(phase_patch(:,:,i))))
            colormap parula;  c = colorbar;
            set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
            set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
            xlabel('\mum'), ylabel('\mum');  axis image
            hold on
            quiver(X,Y,xfield_norm,yfield_norm,'color','w')
            
        end
        
    end

end