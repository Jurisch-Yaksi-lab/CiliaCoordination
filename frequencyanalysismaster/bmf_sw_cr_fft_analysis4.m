function [PowerSpec,Peak,Med,Mean,No_pix,picMask,PeakPhase,mask,nframe, PSD, PeakPos] = bmf_sw_cr_fft_analysis4(data,CBF,varargin)
% Merged Stephan's and Jan's code to both be speedy and generate all plots

% sourceP = Path to the folder where the data is located
% name    = Name of the file containing the data
% Fs      = Frequency of acquisition
% power   = Do you want to retrieve the entire Powerspectrum, true or false.
% targetP = Path to the folder were the data should be saved

% Some information:
% The mask is automatically generated
% Peak frequency is calculated as the max value of the primary frequency map

% Many credits go to Nathalie Jurisch-Yaksi, Jan Niklas Hansen,...
% Stephan Bawolonski, and Benjamin Friedrich for writing this code.

%% Number of frames as small prime factors (not only 2^L)
nframe0=size(data,3); % number of frames should only have small prime factors (2^L best, but not necessary)
nframe=nframe0;
while max(factor(nframe))>5
    nframe=nframe-1;
end

% pixel dimensions of each frame
xmax=size(data,1);
ymax=size(data,2);

% A raw frame for further reference
raw10 = data(:,:,10);

%% Fourier transform ------------------------------------------------------
fprintf('Perform the FFT\n')

% CBF.w_min=15; % [Hz] lower frequency cut-off
iw_min=1+round(CBF.w_min/CBF.Fs*nframe); % corresponding frequency index
nyquist=floor(nframe/2) + 1; % Nyquist frequency OBS: I think nyquist needs a +1

% Allocate memory
PeakPos   = nan(xmax,ymax); % frequency index of peak position
PeakPower = nan(xmax,ymax); % power at peak position
PeakPhase = nan(xmax,ymax); % phase at peak position


PowerSpec = nan(xmax,ymax,nyquist); % entire powerspectrum
tic
for y=1:ymax
    
    % Fourier transform (vectorized for entire pixel column; entire movie would be too large)
    data_fft=squeeze(fft(single(data(:,y,1:nframe)),[],3) );
    
    % Power spectrum
    PowerSpec(:,y,1:nyquist)=abs(data_fft(:,1:nyquist));
    
    % Locate peak in power spectrum: first estimate
    [maxval,maxpos]=max(PowerSpec(:,y,iw_min:end),[],3);
    maxpos=maxpos+iw_min-1;
    
    % Store results
    PeakPos(:,y)   = maxpos;
    PeakPower(:,y) = maxval;
    for x=1:xmax
        PeakPhase(x,y)=angle(data_fft(x,maxpos(x)));
    end
end
toc


%% Retrieving some data

% Translate primary frequency indices (PeakPos) in the power spectra to actual frequencies in Hertz.
%pic = (PeakPos-1).*100./nframe;
pic = ((PeakPos-1).*CBF.Fs)./nframe;

% Create mask
if isempty(varargin)
    picSD=jnh_FreqFilterV1(pic,CBF.SD,CBF.Fs);
    mask = create_maskV2(picSD, CBF.signal_size);
    mask(mask == 0) = NaN;
else
    mask = varargin{1};
end

% Calculate the average powerspectrum
MeanPowerSpec = nanmean(reshape(PowerSpec.*mask,[xmax*ymax,nyquist]));

% Power spectral density 
PSD = 1 /(nframe*CBF.Fs)*PowerSpec.^2;
PSD(2:end-1) = 2*PSD(2:end-1);

% Amplitude
ampl = 2*(PeakPower/ nframe); % Amplitude is just power rescaled/divided by the length of the original time-domain signal

picMask = pic.*mask;
phase = PeakPhase.*mask;
amplMask = ampl.*mask;

% Calculate the Peak, Median, Mean, number of pixel
[~, pos] = max(MeanPowerSpec(iw_min:end));
pos  = pos + iw_min - 1;
Peak = (pos * CBF.Fs) / nframe;
Med = nanmedian(picMask(:));
Mean = nanmean(picMask(:));
No_pix = sum(~isnan(picMask(:)));


% Mask boundary
mask_boundary=~edge(mask,'sobel',[],'nothinning');

% % This is when collecting the peak from the average PowerSpectrum
[ ~, Peak] = max(MeanPowerSpec(iw_min:nyquist));
Peak = ((Peak+iw_min-1)/ nframe)*CBF.Fs;


%% Save data -------------------------------------------------------------
save(fullfile(CBF.targetP,[CBF.name,'_fft.mat']),'PeakPos','PeakPhase','PeakPower', ...
    'mask','nframe','raw10','PowerSpec', 'MeanPowerSpec', 'PSD','-v7.3');

%% Create figures  --------------------------------------------------------

figure; clf
set(gcf,'units','pixels','position',[11 72 1926 1037])

% Raw data frame
ax1=subplot(2,4,1);
imagesc(raw10);
colormap(ax1,'gray'), c=colorbar;
c.Label.String = 'pixel intensity';
title(ax1,'primary freq. peak [mask]')
box off; axis off; axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
% xlabel('\mum'), ylabel('\mum'); 



% Primary frequency peak [raw]
ax2=subplot(2,4,2);
imagesc(pic.*mask_boundary);
colormap(ax2,'jet'), caxis(ax2,CBF.caxis);
c=colorbar; c.Label.String = 'CBF [in Hz]';
title(ax2,'primary freq. peak [raw]')
box off; axis off; axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
% xlabel('\mum'), ylabel('\mum'); 

% Primary frequency peak [mask]
ax3=subplot(2,4,3);
AlphDat =double(~isnan(picMask));
imagesc(picMask, 'AlphaData', AlphDat), colormap(ax3,jet);
caxis(ax3,CBF.caxis), c=colorbar; c.Label.String = 'CBF [in Hz]';
title(ax3,'primary freq. peak [masked]')
box off; axis off; axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
% xlabel('\mum'), ylabel('\mum'); 

% Phase [mask]
ax4=subplot(2,4,4);
AlphDat =double(~isnan(phase));
imagesc(phase, 'AlphaData', AlphDat), colormap(ax4,hsv);
title(ax4,'phase [masked] ');
caxis(ax4,[-pi pi]), c=colorbar; c.Label.String = '\phi (\pi)';
box off; axis off; axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
% xlabel('\mum'), ylabel('\mum'); 

% Freq spectrum mean
ax5=subplot(2,4,5);
plot((2:length(MeanPowerSpec))*(CBF.Fs/nframe),log10(MeanPowerSpec(2:end))*10,'r', 'linewidth', 4) % Messy?
ylabel('Power/frequency (dB/Hz) [mask]')
xlabel ('Frequency (Hz)')
box off; grid off
title(ax5, sprintf('peak frequency = %0.2f', Peak));

% Amplitude --> OBS: make into PSD??
ax6=subplot(2,4,6);
AlphDat =double(~isnan(amplMask));
imagesc(amplMask,'AlphaData', AlphDat)
colormap(ax6,'jet');  colorbar
caxis ([0 max(ampl(:))])
box off, grid off, axis off
title(ax6,'Amplitude');
box off; axis off; axis image
set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
% xlabel('\mum'), ylabel('\mum'); 

% Histogram picMask
ax7=subplot(2,4,7);
histogram(picMask(:))
xlim(CBF.caxis); box off; grid off
title(ax7, sprintf('picMask hist, median = %0.2f, mean = %0.2f', Med, Mean));

% ECDF picMask
ax8=subplot(2,4,8);
try
    ecdf(picMask(:))
catch
    sprintf('There are no values to display for %s', CBF.name)
end
xlim([15 40]); box off; grid off
title(ax8,sprintf('picMask ecdf, no pix = %0.2f',No_pix));

% Put a title without the underscore effect
suptitle(insertBefore(CBF.name,'_','\'));

%% Export nicely
print( '-painters', fullfile(CBF.targetP, CBF.name), '-dpng');
close(gcf)
