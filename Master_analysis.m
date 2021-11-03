%% Master Analysis file
% All analysis performed on light transmission recordings of ciliary beating in the zebrafish nose

%% Step 0: Load the aligned recording 

% Define the path to the raw recording
SourceP = 'C:\Users\chririn\OneDrive - NTNU\Github\03-02-2017_0_0_aligned.mat'; 
CBF.name = '03-02-2017_0_0_aligned'; % Give a name
load(SourceP); 
data = double(aligned); 

% Define the target path 
CBF.targetP = 'C:\Users\chririn\OneDrive - NTNU\Github\03-02-2017_0_0_aligned\';
mkdir(CBF.targetP)
cd(CBF.targetP); 

% Input the frequency of acquisition
CBF.Fs = 107; % Frequency of acquisition

save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
%% Step 1: Defining variables 
 
% Define variables and save them in one common structure
CBF.x = size(data,1); 
CBF.y = size(data,2); 
CBF.w_min = 15; % Lower frequency cutoff 
CBF.caxis = [16 40]; % Upper and lower bound for any frequency plot. 
CBF.spatres = 0.15;  % [um/pixel] spatial resolution

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 2: Perform the fast Fourier transform 

% Fourier Transform (~20s to run the fft & ~200s with plotting)
[PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PSD, PeakPos] = bmf_sw_cr_fft_analysis3(data,CBF);

% It is possible to run the frequency analysis with a predefined mask. 
% [PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PSD, PeakPos] = bmf_sw_cr_fft_analysis3(data,CBF, mask);

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3: Coherence with example pixels 
 
% Decide on reference pixels: uncomment to choose
figure, imagesc(CBF.picSD); 
ref = ginput; 
ref = round(ref);
CBF.ref = ref; % Save the refence pixels

% Set some variables specific to the coherence analysis
CBF.window = hamming(100); % window 100, noverlap 80, nfft 100 work very well! 
CBF.noverlap = 80; 
CBF.nfft = 100; 

% Compute coherence all with one reference pixel (~120s per reference pixel)
[Pxx, val] = cr_coherence_ref(data, CBF);

% Save results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3.5: Coherence versus spectral density scatterplot 
% 
% % Plot coherence versus power spectral contribution (~30s)
% map = cr_coherence_vs_psd(val, PeakPos, PowerSpec, CBF);


%% Step 4: Pairwise Coherence versus distance
% Define the spatial resolution
CBF.f = 5; % subsample factor

% Bin into 100 bins.  
CBF.n = 100;  % Number of bins to segment the powerspectrum into

% May not work outside of Kavlifarm (~ 10min in kavlifarm using the parallel pool). Computationally heavy. 
cr_coherence_all(data, Pxx, CBF)

save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

% %% Step 4.1: Pairwise Coherence versus distance - One patch
% 
% patch_no = 6; 
% % May not work outside of Kavlifarm. Computationally heavy. 
% cr_coherence_all_patch(data, Pxx, CBF, patch_no)
%  

%% Step 5: Frequency segmentation 

% Bin the powerspectrum (<10s)
[freqsBinned] = cr_bin_power_spectrum(CBF, PowerSpec);

% Segment the nose into frequency patches.
CBF.minsize = 400; 
I = freqsBinned.*CBF.mask;

% Run the frequency segmentation (~40s)
[CBF.lmatrix,complist] = identify_frequency_patchesV2(I, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 6: Calculate wavedirection and wavelength. 

% Find the phases for all patches (~50s)
[phase_patch, pos] = calc_phase_patch(data,complist, CBF);

% Calculate the gradient, wave direction and wavelength (~40s)
[CBF.max_wavelength, CBF.bin_size] = cr_patchwise_analysisV2(phase_patch, pos, PowerSpec, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
