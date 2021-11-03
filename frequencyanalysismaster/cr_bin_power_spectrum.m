function [freqsBinned] = cr_bin_power_spectrum(CBF, PowerSpec)
% The resultant powerspectrum has a higher frequency resolution (0.1 Hz)
% than is relevant. 
% This functions groups the powerspectrum into coarser bins and adds the
% powers of the original powerspectrum bins together. 
% CBF = A structure with at least Fs, nframe, n (number of bins), w_min
% (minimum frequency), x, and y. 
% PowerSpec = the powerspectrum per pixel as calculated by the fourier
% transform. 


%% Bin the power spectrum 
freq_res = CBF.Fs / CBF.nframe; % The frequency resolution of the orginal powerspectrum [Hz/frames]
spec = 1:(CBF.nframe/2); 
dCBF = (CBF.nframe/2) / CBF.n; 
iw_min=1+round(CBF.w_min/CBF.Fs*CBF.n); % the position of the minimum frequency

% Define new frequency bins
tp = spec(1):dCBF:spec(end); % All timepoint that will be considered
binEdges = tp(1)-dCBF/2; % Make sure that the bins are centered around the timepoint.
binEdges = [binEdges tp+dCBF/2];

sprintf('The step between frequencies is: %0.2f Hz', dCBF*freq_res)

%% Bin he Power spectrum 

% Allocate for speed
BinSpec = zeros(CBF.x,CBF.y,CBF.n);

% Bin the Power Spectrum accordingly 
for i = 2:CBF.n
    start = floor(binEdges(i)); stop = floor(binEdges(i+1));
    BinSpec(:,:,i) = sum(PowerSpec(:,:,start:stop),3);
end

% Find the peak frequency 
for j = 1:CBF.x
    for k  = 1:CBF.y
        [~, BinFreq(j,k)] = max(squeeze(BinSpec(j,k,iw_min:end)));
    end
end

% Adjust to take the minimum frequency into account. 
BinFreq = BinFreq + iw_min - 1; 

conversion = ((CBF.nframe / 2)/ length(binEdges))*freq_res; % Conversion from old to the new frequency positions 
freqsBinned = BinFreq*conversion; 

%% Plotting figures 

figure, 
subplot(1,2,1), imagesc(CBF.picSD, 'AlphaData', ~isnan(CBF.mask)); c = colorbar; colormap jet, caxis(CBF.caxis)
c.Label.String = 'CBF in [HZ]'; 
title('Original frequency map'); 
axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
xlabel('\mum'), ylabel('\mum'); 

subplot(1,2,2), imagesc(freqsBinned, 'AlphaData', ~isnan(CBF.mask)); c = colorbar; colormap jet, caxis(CBF.caxis)
c.Label.String = 'CBF in [HZ]'; 
title(sprintf('Frequency map: %0.2f Hz binstep', dCBF*freq_res));
axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
xlabel('\mum'), ylabel('\mum'); 

saveas(gcf, [CBF.targetP,CBF.name, sprintf('_figure_binned_freq_%d_bins.png', CBF.n)]);


