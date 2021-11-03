function [PowerSpec,Peak,Med,Mean,No_pix,picMask,PeakPhase,mask,nframe] = bmf_sw_cr_fft_analysis3(data,name,Fs,w_min,targetP)
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
% %% Load data file ---------------------------------------------------------
% fprintf(1,'Load data file %s.\n',name)
%
% % Load matfile
%     mat = matfile(fullfile(sourceP,[name, '.mat']));
%     data = mat.data;
%
% fprintf(1,'Data %s loaded! \n',name)

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

% w_min=15; % [Hz] lower frequency cut-off
iw_min=1+round(w_min/Fs*nframe); % corresponding frequency index
nyquist=floor(nframe/2); % Nyquist frequency

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

% Translate primary frequency indices (PeakPos) in the power spectra to
% actual frequencies in Hertz.
%pic = (PeakPos-1).*100./nframe;
pic = ((PeakPos-1).*Fs)./nframe;

% Create mask
picSD=jnh_FreqFilterV1(pic,3,Fs);
mask = ~isnan(picSD);
mask2 = double(mask); mask2(mask2==0)=NaN;

% Calculate the average powerspectrum
MeanPowerSpec = nanmean(reshape(PowerSpec.*mask2,[xmax*ymax,nyquist]));

% Amplitude
ampl = (PeakPower/ nframe); % Amplitude is just power rescaled/divided by the length of the original time-domain signal

% Doing some masking
mask = double(mask);
mask(mask == 0) = NaN;

picMask = pic.*mask;
phase = PeakPhase.*mask;
amplMask = ampl.*mask;

% Calculate the Peak, Median, Mean, number of pixel
[~, pos] = max(MeanPowerSpec(iw_min:end));
pos  = pos + iw_min - 1;
Peak = (pos * Fs) / nframe;
Med = nanmedian(picMask(:));
Mean = nanmean(picMask(:));
No_pix = sum(~isnan(picMask(:)));


% Mask boundary
mask_boundary=~edge(mask,'sobel',[],'nothinning');

% plot an example PowerSpectrum (just for fun)
% figure, plot(((2:512)/Fs)*nframe, squeeze(PowerSpec(156,348,2:512)/nframe));

% % This is when collecting the peak from the average PowerSpectrum
[ ~, Peak] = max(MeanPowerSpec(iw_min:nyquist));
Peak = ((Peak+iw_min-1)/ nframe)*Fs;


%% Save data -------------------------------------------------------------
save(fullfile(targetP,[name,'fft.mat']),'PeakPos','PeakPhase','PeakPower', ...
    'mask','Fs','nframe', 'name','raw10','PowerSpec', 'MeanPowerSpec','-v7.3');

%% Create figures  --------------------------------------------------------

figure; clf
set(gcf,'units','pixels','position',[11 72 1926 1037])

% Raw data frame
ax1=subplot(2,4,1);
imagesc(raw10);
colormap(ax1,'gray'), c=colorbar;
c.Label.String = 'pixel intensity';
title(ax1,'primary freq. peak [mask]')
box off; axis off

% Primary frequency peak [raw]
ax2=subplot(2,4,2);
imagesc(pic.*mask_boundary);
colormap(ax2,'jet'), caxis(ax2,[15 40]);
c=colorbar; c.Label.String = 'CBF [in Hz]';
title(ax2,'primary freq. peak [raw]')
box off; axis off

% Primary frequency peak [mask]
ax3=subplot(2,4,3);
AlphDat =double(~isnan(picMask));
imagesc(picMask, 'AlphaData', AlphDat), colormap(ax3,jet);
caxis(ax3,[15 40]), c=colorbar; c.Label.String = 'CBF [in Hz]';
title(ax3,'primary freq. peak [masked]')
box off; axis off

% Phase [mask]
ax4=subplot(2,4,4);
AlphDat =double(~isnan(phase));
imagesc(phase, 'AlphaData', AlphDat), colormap(ax4,hsv);
title(ax4,'phase [masked] ');
caxis(ax4,[-pi pi]), c=colorbar; c.Label.String = '\phi (\pi)';
box off; axis off

% Freq spectrum mean
ax5=subplot(2,4,5);
plot((2:length(MeanPowerSpec))*(Fs/nframe),log10(MeanPowerSpec(2:end))*10,'r', 'linewidth', 4) % Messy?
ylabel('Power/frequency (dB/Hz) [mask]')
xlabel ('Frequency (Hz)')
box off; grid off
title(ax5, sprintf('peak frequency = %0.2f', Peak));

% Amplitude
ax6=subplot(2,4,6);
AlphDat =double(~isnan(amplMask));
imagesc(amplMask,'AlphaData', AlphDat)
colormap(ax6,'jet');  colorbar
caxis ([0 max(ampl(:))])
box off, grid off, axis off
title(ax6,'Amplitude');

% Histogram picMask
ax7=subplot(2,4,7);
histogram(picMask(:))
xlim([15 40]); box off; grid off
title(ax7, sprintf('picMask hist, median = %0.2f, mean = %0.2f', Med, Mean));

% ECDF picMask
ax8=subplot(2,4,8);
try
    ecdf(picMask(:))
catch
    sprintf('There are no values to display for %s', name)
end
xlim([15 40]); box off; grid off
title(ax8,sprintf('picMask ecdf, no pix = %0.2f',No_pix));

% Put a title without the underscore effect
suptitle(insertBefore(name,'_','\'));

%% Export nicely
export_fig(fullfile(targetP, name),'-png');
close(gcf)
