function [CBF] = cr_patchwise_analysis(phase_patch, pos, PowerSpec, check, CBF)
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
CBF.NumPatches = size(phase_patch,3); 

% Define this anonymous function 
magnitude = @(x,y) sqrt(x.^2 + y.^2);

close all, 
%% Flip the phase map 
% The y-axis is inverted by default for imagesc() and this creates a discrepancy in the wave direction definition. 
% As a solution we flip the phase and revert the axis on all images. 
% https://stackoverflow.com/questions/19695470/reverse-imagesc-y-axis-tick-marks-to-go-from-positive-to-negative-values

phase_patch = flipud(phase_patch); 
lmatrix = flipud(CBF.lmatrix); 

%% Calculate the wave direction per patch

[xfield, yfield, theta_all, phase_pd, phase_sd, mean_phases, cv_phases, direction_pixel] = wave_direction(lmatrix);  

%% Calculate the wavelength per patch

[wave_length, magn_all] = wavelength(lmatrix);

save([CBF.targetP, CBF.name, '_results_wave_direction_and_lengths.mat'], 'wave_length', 'phase_pd', 'phase_sd', 'direction_pixel','-v7.3');
close all, 

%% Run some checks for every patch 

if check == true
    cr_check_per_patch(lmatrix, phase_patch, xfield, yfield, phase_pd, mean_phases, cv_phases, pos, PowerSpec, theta_all, magn_all, wave_length, CBF)
else
    sprintf('No checks were run')
end

%% Nested functions
    function [xfield, yfield, theta_all, phase_pd, phase_sd, mean_phases, cv_phases, direction_pixel] = wave_direction(lmatrix)
     
        % Allocate for speed
        mean_phases = nan(1,CBF.NumPatches); % Circular mean per patch
        cv_phases  = nan(1,CBF.NumPatches); % Circular variance per patch
        phase_pd = nan(size(lmatrix)); % To store the new patch image that encodes mean directions
        phase_sd = nan(size(lmatrix)); % To store the new patch image that encodes the variance of directions
        direction_pixel = nan(size(lmatrix));  % To store pixel-wise directions
        xfield = nan([size(lmatrix),CBF.NumPatches]);
        yfield = nan([size(lmatrix),CBF.NumPatches]);
        theta_all = cell(1,CBF.NumPatches); 
        
        %% Determine gradient directions for the frequency patch
        
        for iPatch = 1:CBF.NumPatches
            
            % Determine directions of gradient vectors
            mask_patch = lmatrix == iPatch;
            [theta_patch, xfield(:,:,iPatch), yfield(:,:,iPatch)] = determine_theta(phase_patch(:,:,iPatch), 1);
          
            % Save the shape of theta_patch
            dirs = wrapTo2Pi(theta_patch); 
            
            % Restrict directions to the frequency patch we want to investigate
            theta_patch = theta_patch(~isnan(theta_patch)); 
            
            % Determine mean phase angle
            if ~isempty(theta_patch)
                mean_phases(iPatch) = angle(sum(exp(1i.*theta_patch(:)))); 
                cv_phases(iPatch) = 1-abs(mean(exp(1i.*theta_patch(:)))); % determine circular standard deviation of mean phase angle
            else
                mean_phases(iPatch) = nan;
                cv_phases(iPatch) = nan;
            end
            
            % Store the mean direction in a frequency patch image
            phase_pd(mask_patch) = wrapTo2Pi(mean_phases(iPatch)); % Wrap from [0 to 2pi]
            phase_sd(mask_patch) = cv_phases(iPatch);
            direction_pixel(mask_patch) = dirs(mask_patch);
            theta_all{iPatch} = dirs; 
            theta_patch = []; %#ok<NASGU> % Clear the variable to avoid overlap
                                  
        end

        %% Create figures
        
%         % Phase direction map
%         figure,
%         imagesc(phase_pd, 'AlphaData', ~isnan(phase_pd))
%         colormap hsv, hold on
%         set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
%         set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
%         set(gca,'YDir','normal');
%         xlabel('\mum'), ylabel('\mum'); axis image,  caxis([0 2*pi])
%         saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_direction.png']);
%         % print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_direction'], '-depsc');
        
%         % Phase direction variance map
%         figure,
%         imagesc(phase_sd, 'AlphaData', ~isnan(phase_pd))
%         colormap gray,  c = colorbar;
%         hold on, caxis([0 1]), c.Label.String = 'standard deviation';
%         set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
%         set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
%         set(gca,'YDir','normal');
%         xlabel('\mum'), ylabel('\mum'); axis image
%         saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_direction_variance.png']);
        
        % Phase per patch map
        figure,
        imagesc(sum(phase_patch,3, 'omitnan'), 'AlphaData', ~isnan(phase_pd))
        colormap hsv,  c = colorbar; caxis([0 2*pi])
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        set(gca,'YDir','normal');
        xlabel('\mum'), ylabel('\mum'); axis image
        hold on
        c.Label.String = 'phase angle';
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_phasemap.png']);
        print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_phasemap'], '-depsc');
        
        % Phase direction and variance map
        figure,
        imagesc(phase_pd, 'AlphaData', 1 - phase_sd)
        set(gca,'YDir','normal');
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        colormap hsv,  c = colorbar;  caxis([0 2*pi])
        c.Label.String = 'wave direction'; axis image
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_patch_combined_direction_and_variance.png']);
        print( '-painters', [CBF.targetP, CBF.name, '_figure_patch_combined_direction_and_variance'], '-depsc');   
        
    end

    function [wave_length, magn_all] = wavelength(lmatrix)
                
        % Allocate for speed
        max_lambdas = nan(1, CBF.NumPatches);
        magn_rad = nan(size(lmatrix)); 
        magn_pix = nan(size(lmatrix)); 
        wave_length = nan(size(lmatrix));
        magn_all = cell(1,CBF.NumPatches); 
        
        %% Determine gradient magnitudes for the frequency patch
        
        for iPatch = 1:CBF.NumPatches
            
            % Determine magnitudes
            magn = magnitude(squeeze(xfield(:,:,iPatch)), squeeze(yfield(:,:,iPatch))); % [rad/pixel]
            
            % Determine metachronal wave lengths
            magn_patch = magn(~isnan(magn));
            lambdas = 2*pi./magn_patch; % Wavelength in [pixels] for all pixels
            mask_patch = ~isnan(magn);
            
            % Store the mean direction in a frequency patch image
            magn_rad(mask_patch) = magn_patch; % [rad/pixel]
            magn_pix(mask_patch) = lambdas; % [pixel]
            
            % Attempt to identify one common wavelength per patch (sometimes a patch is empty)
            try
                % Plot the histogram of lambdas and overlay a kernel
                h = histfit(lambdas(lambdas<CBF.max_wavelength),CBF.bin_size, 'kernel');
                
                % Find the peak of the kernel by taking the max value of the kernel
                [~,p] = max(h(2).YData);
                max_lambdas(iPatch) = h(2).XData(p); % Extract the value at the peak.

                % Save the wavelength
                wave_length(lmatrix==iPatch) = max_lambdas(iPatch)*CBF.spatres; % wavelength in [µm] per patch       
            catch
                sprintf('No lambdas found in patch %d', iPatch);
            end
                    
            % Save for plotting later
            magn_all{iPatch} = magn_rad; 
            magn_rad = nan(size(lmatrix));
        end
        
        %% Create figures
        
%         % Magnitude in [rad/pixels]
%         figure, imagesc(magn_rad, 'AlphaData', ~isnan(magn_rad));
%         set(gca,'YDir','normal');
%         c = colorbar;
%         c.Label.String = 'Magnitude in [rad/pixel]';
%         set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
%         set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
%         xlabel('\mum'), ylabel('\mum')
%         axis image
%         hold on
%         caxis([0 1.8])
%         title('Image gradient magnitude per patch');
%         saveas(gcf, [CBF.targetP, CBF.name, '_figure_magnitude_rad_pixel.png']);
        
%         % Magnitude in [pixels]
%         figure, imagesc(magn_pix, 'AlphaData', ~isnan(magn_pix));
%         set(gca,'YDir','normal');
%         c = colorbar;
%         c.Label.String = 'Magnitude in [pixel]';
%         set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
%         set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
%         xlabel('\mum'), ylabel('\mum');  axis image
%         hold on
%         caxis([0 500])
%         title('Image gradient magnitude per patch');
%         saveas(gcf, [CBF.targetP, CBF.name, '_figure_magnitude_pixel.png']);
        
        % Wavelength in micrometer per patch
        figure, imagesc(wave_length, 'AlphaData', ~(lmatrix==0));
        set(gca,'YDir','normal');
        c = colorbar; colormap jet;
        c.Label.String = 'Wavelength [µm]';
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum');  axis image
        hold on
        caxis([0 8])
        title('wavelength per patch');
        saveas(gcf, [CBF.targetP, CBF.name, '_figure_wavelength_per_patch.png']);
        print( '-painters', [CBF.targetP, CBF.name, '_figure_wavelength'], '-depsc');
               
    end

end