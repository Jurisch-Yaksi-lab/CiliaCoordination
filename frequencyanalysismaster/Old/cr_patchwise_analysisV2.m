function [] = cr_patchwise_analysisV2(w_min)
% In this function, a gradient map is calculated for a static phase. 
% From this gradient map, both the wave direction and wavelength will be
% calculated. 

% TODO: 
% which phase map will I use for the calculation. Which one I use will
% probably matter... 


%% Variables 

% w_min = 15; % in Hz
FreqRes = Fs / nframe; 
iw_min = round(w_min/FreqRes);
spatres = 0.15;  % [um/pixel] spatial resolution
minsize = 400; 
%% Script for wave direction per patch
magnitude = @(x,y) sqrt(x.^2 + y.^2);

% double_mean = @(I) squeeze(nanmean(nanmean(I))); % May come in handy
% %% One frequency per patch --> one calculated phase 
% 
% % Find the phases for all patches 
% [phase_patch, pos] = calc_phase_patch(data,Fs,lmatrix,complist);
% phase_patch = phase_patch+pi; % shift phases to interval [0, 2*pi]

NumPatches = size(phase_patch,3); 

%%  Double-check whether the frequency of the phase-calculation fits the frequency calculated elsewhere. 
% 
% figure(1), 
% for n = 1:NumPatches
%     % non-binned frequency --> I don't know how binning will affect the phases
%     patch_ps = double_mean(apply_mask(lmatrix == n,PowerSpec));
%     subplot(5,3,n), plot(patch_ps(3:end)), ylim([0 500]); title(sprintf('patch, %i', n));
%     [~, id] = max(patch_ps(iw_min:end));
%     pos2(n) = (id+iw_min-1); % Makes sense :)
% end

%% Calculate the wave direction with these new wave-per-patch calculations
mean_phases = nan(1,NumPatches);
cv_phases  = nan(1,NumPatches);
binEdges = 0:(2*pi)/12:2*pi;
phase_pd = nan(size(lmatrix)); % to store the new patch image that encodes mean directions
phase_sd = nan(size(lmatrix)); 
 


for i = 1:NumPatches
          
    % determine directions of gradient vectors
    mask2 = lmatrix == i; 
    [theta_patch, xfield(:,:,i), yfield(:,:,i)] = determine_theta(phase_patch(:,:,i), 1);
    
    % restrict directions to the frequency patch we want to investigate
    theta_patch = theta_patch(~isnan(theta_patch));
    
    % determine mean phase angle
    if ~isempty(theta_patch)
        mean_phases(i) = angle(sum(exp(1i.*theta_patch(:))));
    % determine circular standard deviation of mean phase angle
         cv_phases(i) = 1-abs(mean(exp(1i.*theta_patch(:))));
    else
        mean_phases(i) = nan;
        cv_phases(i) = nan;
    end
     
    % encode the mean direction in a frequency patch image
    phase_pd(mask2) = mean_phases(i)+pi; % shift the phases to [0, 2pi] interval
    phase_sd(mask2) = cv_phases(i); 
    
    figure(2), subplot(6,5,i), polarhistogram(theta_patch+pi,binEdges), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));% xlim([0,2*pi])
    figure(3), subplot(6,5,i), imagesc(phase_pd, 'AlphaData', mask2), title(sprintf('mean_phase = %0.2f, sv = %0.2f', mean_phases(i)+pi, cv_phases(i)));% xlim([0,2*pi]), 
    colormap hsv
end

%% create figure
figure(4),
imagesc(phase_pd, 'AlphaData', ~isnan(phase_pd))
colormap hsv
hold on

set(gca,'XTickLabel',get(gca,'XTick')*spatres)
set(gca,'YTickLabel',get(gca,'YTick')*spatres)
xlabel('\mum')
ylabel('\mum')

set(gcf,'color','w')
set(gcf,'units','pixels','position',[679  587  560  420]) 
    
figure(5),
imagesc(phase_sd, 'AlphaData', ~isnan(phase_pd))
colormap gray,  c = colorbar; 
hold on
caxis([0 1]),
c.Label.String = 'standard deviation';

set(gca,'XTickLabel',get(gca,'XTick')*spatres)
set(gca,'YTickLabel',get(gca,'YTick')*spatres)
xlabel('\mum')
ylabel('\mum')

set(gcf,'color','w')
set(gcf,'units','pixels','position',[679  587  560  420]) 

figure(6),
imagesc(nansum(phase_patch,3), 'AlphaData', ~isnan(phase_pd))
colormap hsv,  c = colorbar; 
hold on
c.Label.String = 'phase angle';

set(gcf,'units','pixels','position',[679  587  560  420]) 


figure(10),
imagesc(phase_pd, 'AlphaData', 1 - phase_sd)
colormap hsv,%  c = colorbar; 
% hold on
% c.Label.String = 'wave direction';
set(gcf,'units','pixels','position',[679  587  560  420]) 


%% Plot the heatmap for that particular frequency 
figure(7),
for n = 1:NumPatches
map(:,:,n) = squeeze(PowerSpec(:,:,pos(n))); 

% Plot the heatmap 
borders = edge(lmatrix==n, 'canny');
subplot(4,4,n), imagesc(squeeze(map(:,:,n)), 'AlphaData', ~borders), colormap hot; c = colorbar, hold on, 
title(sprintf('heatmap for patch %i', n))
c.Label.String = 'spectral power';
caxis([0 1800])
end

%% Calculate wave length 

% Shift phases to interval [0, 2*pi] since otherwise our filter won't work
% properly.
 % shift phases to interval [0, 2*pi]

%% preprocess data
% determine areas for all frequency clusters
stats = regionprops(complist,'Area');
% areas = [stats.Area]*spatres^2;
areas = [stats.Area]; % Describe it in pixels for now

%% Determine gradient magnitudes for the frequency patch

% % Work through all the frequency patches of minimum size
% numPatches = sum(areas>minsize);

% allocate memory
median_lambda = nan(1, NumPatches);
% patchIdx = find(areas>minsize);
lambda_all = [];

for i = 1:NumPatches
    % Determine magnitudes
    magn = magnitude(squeeze(xfield(:,:,i)), squeeze(yfield(:,:,i)));
    
    % restrict to a single frequency patch
    mask2 = lmatrix == i;
    magn_patch = mask_image(magn,mask2);
    m = magn_patch;
    
    % determine metachronal wave lengths
    magn_patch = magn_patch(~isnan(magn_patch));
    lambdas = 2*pi./magn_patch;
    lambda_all = [ lambda_all lambdas'];
    % determine median and convert it to micrometers
    median_lambda(i) = median(lambdas(:));%*spatres;
    
%     figure, hist(lambdas(lambdas<=100)), hold on, 
%     figure, hist(lambdas), hold on, 
%     plot([median_lambda(i),median_lambda(i)], [0 1400])
%      saveas(gcf, ['X:\Christa\ANALYZED\Data\19b\Tasks for 27.02.2020\Patch_wise analysis\WaveLength\', sprintf('HistOfWavelengths2 patch%i.png', i)]);
    figure, imagesc(m, 'AlphaData', mask2),  colormap jet,title('Magnitude per patch'), colorbar%, caxis([0 5])
      saveas(gcf, ['X:\Christa\ANALYZED\Data\19b\Tasks for 27.02.2020\Patch_wise analysis\WaveLength\', sprintf('Magnitude patch%i.png', i)]);
 close all, 
end

% close all
%% Plot the figure

mlambdas = nan(size(lmatrix)); % allocate memory
for i = 1:NumPatches
    mask2 = lmatrix == i;
    mlambdas(mask2) = median_lambda(i);
end

mask3 = ~isnan(mlambdas); 
figure(8)
imagesc(mlambdas*spatres, 'AlphaData', mask3)
colormap jet
c = colorbar;
c.Label.String = 'WaveLength in µm';
hold on
caxis([0 18])

%% Determine gradient magnitudes for the frequency patch --> Now adjusting the patch size

% % Work through all the frequency patches of minimum size
% numPatches = sum(areas>minsize);

% allocate memory
median_lambda = nan(1, NumPatches);
% patchIdx = find(areas>minsize);
lambda_all = [];

for i = 1:NumPatches
    
   % restrict to a single frequency patch
   % and imdilate to get rid of the small border-vectors
    
   
       % Determine magnitudes
%      magn = magnitude(squeeze(xfield(:,:,i)), squeeze(yfield(:,:,i)));
    magn = squeeze(xfield(:,:,i))+ squeeze(yfield(:,:,i));

   mask2 = lmatrix == i;
    magn_patch = mask_image(magn,mask2);
%     m = magn_patch;
%     mask3 = imerode(mask2,strel('disk', 2)); % This kind of works
%     figure, imagesc(m, 'AlphaData', mask3)
%     figure, imagesc(m, 'AlphaData', (~mask3 & mask2) )
    
    
     % determine metachronal wave lengths
    magn_patch = magn_patch(~isnan(magn_patch));
    lambdas = (2*pi)./magn_patch;
    
    lambda_all = [ lambda_all lambdas'];
    % determine median and convert it to micrometers
    median_lambda(i) = nanmedian(lambdas(:));%*spatres;
    
%     figure, hist(lambdas(lambdas<=100)), hold on, 
%     lambdas_plot = lambdas(lambdas<200);
    figure, histogram(lambdas,70), hold on, 
    figure, imagesc((2*pi)./magn, 'AlphaData', (mask2&~isnan(magn)))
%     plot([median_lambda(i),median_lambda(i)], [0 1400])
%      saveas(gcf, ['X:\Christa\ANALYZED\Data\19b\Tasks for 27.02.2020\Patch_wise analysis\WaveLength\', sprintf('HistOfWavelengths2 patch%i.png', i)]);
%     figure, imagesc(m, 'AlphaData', mask2),  colormap jet,title('Magnitude per patch'), colorbar%, caxis([0 5])
%       saveas(gcf, ['X:\Christa\ANALYZED\Data\19b\Tasks for 27.02.2020\Patch_wise analysis\WaveLength\', sprintf('Magnitude patch%i.png', i)]);
%  close all, 
end

% close all
%% Plot the figure

mlambdas = nan(size(lmatrix)); % allocate memory
for i = 1:NumPatches
    mask2 = lmatrix == i;
    mlambdas(mask2) = median_lambda(i);
end

mask3 = ~isnan(mlambdas); 
figure(8)
imagesc(mlambdas.*spatres, 'AlphaData', mask3)
colormap jet
c = colorbar;
c.Label.String = 'WaveLength in µm';
hold on
